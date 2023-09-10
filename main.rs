use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::sync::{Arc, RwLock};
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use clap::{Arg, App};
use csv;

#[derive(Debug, Clone)]
struct GeneInfo {
    contig_id: String,
    start_position: usize,
    end_position: usize,
    gene_id: String,
}

#[derive(Debug)]
enum GeneReplaceError {
    CsvError(csv::Error),
    IoError(io::Error),
}

impl std::error::Error for GeneReplaceError {}

impl std::fmt::Display for GeneReplaceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GeneReplaceError::CsvError(err) => write!(f, "CSV error: {}", err),
            GeneReplaceError::IoError(err) => write!(f, "IO error: {}", err),
        }
    }
}

impl From<GeneReplaceError> for io::Error {
    fn from(error: GeneReplaceError) -> Self {
        match error {
            GeneReplaceError::CsvError(csv_err) => io::Error::new(io::ErrorKind::Other, format!("CSV error: {}", csv_err)),
            GeneReplaceError::IoError(io_err) => io_err,
        }
    }
}

fn gene_snv_replace(
    contigs_file: &str,
    mutations_file: &str,
    gene_positions_file: &str,
    output_file: &str,
    _gene_contigs_file: &str,
    num_threads: usize,
) -> io::Result<HashMap<String, String>> {
    // 读取 contigs、mutations 和 gene positions
    let contigs = read_contigs(contigs_file)?;
    let mutations = read_mutations(mutations_file)?;

    // 读取基因位置信息
    let gene_positions_map = read_gene_positions(gene_positions_file).map_err(|e| {
        eprintln!("Error reading gene positions: {:?}", e);
        io::Error::new(io::ErrorKind::Other, "Gene position reading error")
    })?;

    // 创建线程池
    let pool = ThreadPoolBuilder::new().num_threads(num_threads).build().unwrap();

    // 使用 Arc 和 RwLock 创建存储处理结果的 HashMap
    let mutated_genes: Arc<RwLock<HashMap<String, String>>> = Arc::new(RwLock::new(HashMap::new()));

    // 在并行处理突变之前，先组织基因信息，以基因ID为键，对应基因信息及其对应的contigs列表为值
    let mut gene_contigs_map: HashMap<String, Vec<GeneInfo>> = HashMap::new();
    for gene_info in gene_positions_map.values().flatten() {
        gene_contigs_map.entry(gene_info.gene_id.clone()).or_insert(vec![]).push(gene_info.clone());
    }

    // 并行处理突变
    pool.install(|| {
        let mutated_genes_clone = Arc::clone(&mutated_genes);
        gene_contigs_map.par_iter().for_each(|(gene_id, gene_info_list)| {
            for gene_info in gene_info_list {
                let contig_id = &gene_info.contig_id;
                if let Some(contig_sequence) = contigs.get(contig_id) {
                    let contig_sequence = contig_sequence.to_owned();
                    let mutated_sequence = {
                        let mut mutated_contig = contig_sequence;
                        for position in gene_info.start_position..=gene_info.end_position {
                            if let Some(new_base) = mutations.iter().find(|(contig, pos, _)| contig == contig_id && *pos == position).map(|(_, _, base)| *base) {
                                let gene_position = position - 1;
                                if gene_position < mutated_contig.len() {
                                    let mut chars: Vec<char> = mutated_contig.chars().collect();
                                    chars[gene_position] = new_base;
                                    mutated_contig = chars.iter().collect();
                                }
                            }
                        }
                        mutated_contig
                    };
                    let gene_id_clone = gene_id.clone();
                    (*mutated_genes_clone.write().unwrap()).insert(gene_id_clone, mutated_sequence.clone());
                } else {
                    eprintln!("找不到contigs序列：{}", contig_id);
                }
            }
        });
    });

    // 接下来我们需要根据基因id和contigs_id的对应关系，从mutated_genes中提取对应的突变后基因序列
    let mutated_genes_read = mutated_genes.read().unwrap();
    let mutated_genes_result: HashMap<String, String> = gene_contigs_map
        .iter()
        .filter_map(|(gene_id, gene_info_list)| {
            if let Some(mutated_sequence) = mutated_genes_read.get(gene_id) {
                // 获取基因对应的所有contigs信息列表
                for gene_info in gene_info_list {
                    let contig_id = &gene_info.contig_id;
                    if let Some(_contig_sequence) = contigs.get(contig_id) {
                        // 截取基因序列的起始和终止位置
                        let start = gene_info.start_position - 1; // 起始位点（因为序列索引从0开始）
                        let end = gene_info.end_position; // 终止位点
                        let gene_sequence = &mutated_sequence[start..end];
                        // 将基因ID和突变后的基因序列存储在哈希表中
                        return Some((gene_id.clone(), gene_sequence.to_string()));
                    } else {
                        eprintln!("找不到contigs序列：{}", contig_id);
                    }
                }
            } else {
                eprintln!("找不到对应突变基因：{}", gene_id);
            }
            None
        })
        .collect();

    // 在这里我们修改代码，将突变后的基因序列写入到文件中（创建新文件）
    let mut output_file = File::create(output_file)?;
    for (gene_id, mutated_sequence) in mutated_genes_result.iter() {
        writeln!(output_file, ">{}\n{}", gene_id, mutated_sequence)
            .map_err(|e| {
                eprintln!("Error writing to output file: {}", e);
                io::Error::new(io::ErrorKind::Other, "Output file writing error")
            })?;
    }

    Ok(mutated_genes_result)
}


fn main() -> io::Result<()> {
    let matches = App::new("z10")
        .arg(Arg::with_name("contigs_file")
            .required(true)
            .takes_value(true)
            .index(1)
            .help("Contigs 文件路径"))
        .arg(Arg::with_name("mutations_file")
            .required(true)
            .takes_value(true)
            .index(2)
            .help("突变信息文件路径"))
        .arg(Arg::with_name("gene_positions_file")
            .required(true)
            .takes_value(true)
            .index(3)
            .help("基因位置文件路径"))
        .arg(Arg::with_name("output_file")
            .required(false)
            .takes_value(true)
            .short("o")
            .long("output")
            .help("输出文件路径"))
        .arg(Arg::with_name("gene_contigs_file")
            .required(false)
            .takes_value(true)
            .short("gc")
            .long("gene-contigs")
            .help("基因和contigs对应关系文件路径"))
        .arg(Arg::with_name("num_threads")
            .required(false)
            .takes_value(true)
            .short("t")
            .long("num_threads")
            .help("线程数"))
        .get_matches();

    let contigs_file = matches.value_of("contigs_file").unwrap_or_else(|| {
        eprintln!("未提供 Contigs 文件路径！");
        std::process::exit(1);
    });

    let mutations_file = matches.value_of("mutations_file").unwrap_or_else(|| {
        eprintln!("未提供突变信息文件路径！");
        std::process::exit(1);
    });

    let gene_positions_file = matches.value_of("gene_positions_file").unwrap_or_else(|| {
        eprintln!("未提供基因位置文件路径！");
        std::process::exit(1);
    });
    
    let output_file = matches.value_of("output_file").unwrap_or("output.fasta");
    println!("输出突变基因序列: {}", output_file);
    let gene_contigs_file = matches.value_of("gene_contigs_file").unwrap_or("gene_contigs.txt");
    let num_threads: usize = matches
        .value_of("num_threads")
        .and_then(|val| val.parse().ok())
        .unwrap_or_else(num_cpus::get);

    // 调用 gene_snv_replace 函数并获取 mutated_genes 的结果
    let mutated_genes = gene_snv_replace(
        contigs_file,
        mutations_file,
        gene_positions_file,
        output_file,
        gene_contigs_file,
        num_threads,
    )?;

    // 在这里可以处理 mutated_genes 变量，比如输出到控制台或保存到文件等
    for (_gene_id, _mutated_sequence) in mutated_genes.iter() {

    }
    let _current_dir = std::env::current_dir().unwrap();


    // 手动刷新 stdout，确保立即显示输出
    std::io::stdout().flush().unwrap();

    Ok(())
}

fn read_contigs(filename: &str) -> io::Result<HashMap<String, String>> {

    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut contigs = HashMap::new();
    let mut current_id = String::new();
    let mut current_sequence = String::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            if !current_id.is_empty() && !current_sequence.is_empty() {
                contigs.insert(current_id.clone(), current_sequence.clone());
            }
            current_id = line[1..].to_string();
            current_sequence.clear();
        } else {
            current_sequence.push_str(&line);
        }
    }

    if !current_id.is_empty() && !current_sequence.is_empty() {
        contigs.insert(current_id.clone(), current_sequence.clone());
    }

    Ok(contigs)
}

fn read_mutations(filename: &str) -> io::Result<Vec<(String, usize, char)>> {

    let file = File::open(filename)?;
    let reader = BufReader::new(file);

    let mut mutations = vec![];

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() == 3 {
            let contig_id = parts[0].to_string();
            if let Ok(position) = parts[1].parse::<usize>() { // 位置信息在第3个部分
                if let Some(new_base) = parts[2].chars().next() { // 新碱基在第4个部分
                    mutations.push((contig_id, position, new_base));
                }
            }
        }
    }
    // 在函数末尾添加调试输出
    std::io::stdout().flush().unwrap(); // 刷新输出缓冲区
    Ok(mutations)
}

fn read_gene_positions(filename: &str) -> io::Result<HashMap<String, Vec<GeneInfo>>> {
    let file = std::fs::File::open(filename).map_err(GeneReplaceError::IoError)?;
    let mut rdr = csv::ReaderBuilder::new().has_headers(false).delimiter(b',').from_reader(file);

    let mut gene_positions_map: HashMap<String, Vec<GeneInfo>> = HashMap::new();

    for result in rdr.records() {
        let record = result.map_err(GeneReplaceError::CsvError)?;
        let record_data = record.iter().map(|field| field.trim()).collect::<Vec<_>>();
        if record_data.len() == 4 {
            let gene_id = record_data[1].to_string();
            let contig_id = record_data[0].to_string();
            if let Ok(start_position) = record_data[2].parse::<usize>() {
                if let Ok(end_position) = record_data[3].parse::<usize>() {
                    let gene_info = GeneInfo {
                        contig_id: contig_id.clone(),
                        start_position,
                        end_position,
                        gene_id: gene_id.clone(),
                    };
                    gene_positions_map.entry(contig_id).or_insert(vec![]).push(gene_info);
                } else {
                    eprintln!("Error parsing end_position field");
                }
            } else {
                eprintln!("Error parsing start_position field");
            }
        } else {
            eprintln!("Invalid row format: {:?}", record);
        }
    }
    Ok(gene_positions_map)
}
