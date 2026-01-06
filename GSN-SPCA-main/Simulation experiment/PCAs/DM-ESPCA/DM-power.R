setwd("I:\\ESPCA\\模拟实验\\定稿\\DM-ESPCA\\DM-ESPCA")
library(readr)
library(dplyr)

process_gene_data <- function(input_file, output_file) {
  # 增强型文件读取（处理最后一行无换行符问题）
  raw_data <- readLines(input_file, warn = FALSE) %>% 
    trimws() %>% 
    .[nzchar(.)]  # 过滤空行
  
  # 初始化结果向量
  results <- vector("numeric", length(raw_data))
  
  # 添加进度条
  pb <- txtProgressBar(min = 0, max = length(raw_data), style = 3)
  
  for (i in seq_along(raw_data)) {
    line <- raw_data[i]
    
    # 数据清洗（带有效性检查）
    data <- suppressWarnings(
      as.numeric(unlist(strsplit(line, "\\s+")))
    ) %>% 
      na.omit() %>% 
      .[is.finite(.)]  # 过滤无穷值
    
    # 有效性检查（生物数据合理性）
    if (length(data) < 2) {
      results[i] <- 1.5
      next
    }
    
    # 动态分组策略（防止奇数长度问题）
    split_point <- max(1, floor(length(data)/2))
    group1 <- data[1:split_point]
    group2 <- data[(split_point+1):length(data)]
    
    # 稳健统计检验（处理全零方差情况）
    p_value <- tryCatch({
      if (sd(group1) == 0 && sd(group2) == 0) {
        1.0  # 两组均为常数
      } else {
        test <- t.test(group1, group2, var.equal = FALSE)
        pmax(pmin(test$p.value, 1.0), 0.0)
      }
    }, error = function(e) 1.0)
    
    # 边界安全的归一化计算
    results[i] <- round(2 - p_value, 4) %>% 
      pmax(1.0) %>% 
      pmin(2.0)
    
    # 更新进度
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # 写入结果文件（兼容科学计数法）
  write_lines(format(results, nsmall = 4, scientific = FALSE), output_file)
  
  # 返回处理统计信息
  list(
    total_lines = length(raw_data),
    valid_lines = sum(results != 1.5),
    success_rate = mean(results != 1.5)
  )
}

# 执行处理并获取元数据
meta_info <- process_gene_data("gene_new_o.txt", "power.txt")

# 打印处理摘要
cat(sprintf("\n处理完成：\n总行数：%d\n有效行：%d\n成功率：%.1f%%\n",
            meta_info$total_lines,
            meta_info$valid_lines,
            meta_info$success_rate * 100))