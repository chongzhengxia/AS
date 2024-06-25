rm(list = ls())
if(TRUE){
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#读取A3SS合并后的文件，添加identity列
A3SS <- read.table(file = "C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\head_all_A3SS.MATS.JC.txt",header = TRUE)
A3SS$Identity <- NA
for (number_row in 1:nrow(A3SS)){
  A3SS[number_row,25] <- paste(as.character(A3SS[number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
}
#unique事件
Identity <- unique(sort(A3SS$Identity))
A3SS_psi_matrix <- data.frame(Identity = Identity)
#循环读取文件
filenames <- c("H9DW0030","H9DW0038","H9DW0039","H9DW0047","H9DW0049","H9DW0052","H9DW0056","H9DW0057","H9DW0064","H9DW0068","H9DW0073","H9DW0074","H9DW0076","H9DW0085","H9DW0092","H9DW0109","H9DW0110","H9DW0116","H9DW0117","H9DW0123","H9DW0125","H9DW0126","H9DW0154","H9DW0169","H9DW0192","H9DW0193","H9DW0194","H9DW0201","H9DW0212","H9DW0219","H9DW0227","H9DW0230","H9DW0250","H9DW0251","H9DW0261","H9DW0262","H9DW0273","H9DW0274","H9DW0280","H9DW0300","H9DW0301","H9DW0312","H9DW0314","H9DW0330","H9DW0338","H9DW0344","H9DW0348","H9DW0351","H9DW0362","H9DW0368","H9DW0377","H9DW0378","H9DW0380","H9DW0383","H9DW0385","H9DW0388","H9DW0079","H9DW0327")
dirpath <- 'C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\A3SS'
result_A3SS <- list()
for (filename in filenames){
  df <- read.table(paste(dirpath,paste(filename,'_A3SS.MATS.JC.txt',sep = ''),sep = "\\"),header = TRUE)
  result_A3SS[[paste(filename,'_A3SS.MATS.JC.txt',sep = '')]] <- df
}
#将每个文件的inclevel改名为psi并且添加identity列
for(i in names(result_A3SS)){
  colnames(result_A3SS[[i]])[21] <- paste(paste(substr(i,1,8),'N',sep = ''),'_PSI',sep = '')
  colnames(result_A3SS[[i]])[22] <- paste(paste(substr(i,1,8),'T',sep = ''),'_PSI',sep = '')
  result_A3SS[[i]]$Identity <- NA
  for (number_row in 1:nrow(result_A3SS[[i]])){
    result_A3SS[[i]][number_row,25] <- paste(as.character(result_A3SS[[i]][number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
  }
}
#取出文件的psi的两列以及Identity列，进行循环合并
for(i in names(result_A3SS)){
  tmp_only_psi_dataframe <- result_A3SS[[i]][,c(25,21,22)]
  A3SS_psi_matrix <- merge(A3SS_psi_matrix,tmp_only_psi_dataframe,by = "Identity",all.x  = TRUE)
}
write.table(A3SS_psi_matrix, "C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\A3SS_psi_matrix.txt", sep="\t",quote = FALSE, row.names = FALSE)
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#读取A5SS合并后的文件，添加identity列
A5SS <- read.table(file = "C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\head_all_A5SS.MATS.JC.txt",header = TRUE)
A5SS$Identity <- NA
for (number_row in 1:nrow(A5SS)){
  A5SS[number_row,25] <- paste(as.character(A5SS[number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
}
#unique事件
Identity <- unique(sort(A5SS$Identity))
A5SS_psi_matrix <- data.frame(Identity = Identity)
#循环读取文件
filenames <- c("H9DW0030","H9DW0038","H9DW0039","H9DW0047","H9DW0049","H9DW0052","H9DW0056","H9DW0057","H9DW0064","H9DW0068","H9DW0073","H9DW0074","H9DW0076","H9DW0085","H9DW0092","H9DW0109","H9DW0110","H9DW0116","H9DW0117","H9DW0123","H9DW0125","H9DW0126","H9DW0154","H9DW0169","H9DW0192","H9DW0193","H9DW0194","H9DW0201","H9DW0212","H9DW0219","H9DW0227","H9DW0230","H9DW0250","H9DW0251","H9DW0261","H9DW0262","H9DW0273","H9DW0274","H9DW0280","H9DW0300","H9DW0301","H9DW0312","H9DW0314","H9DW0330","H9DW0338","H9DW0344","H9DW0348","H9DW0351","H9DW0362","H9DW0368","H9DW0377","H9DW0378","H9DW0380","H9DW0383","H9DW0385","H9DW0388","H9DW0079","H9DW0327")
dirpath <- 'C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\A5SS'
result_A5SS <- list()
for (filename in filenames){
  df <- read.table(paste(dirpath,paste(filename,'_A5SS.MATS.JC.txt',sep = ''),sep = "\\"),header = TRUE)
  result_A5SS[[paste(filename,'_A5SS.MATS.JC.txt',sep = '')]] <- df
}
#将每个文件的inclevel改名为psi并且添加identity列
for(i in names(result_A5SS)){
  colnames(result_A5SS[[i]])[21] <- paste(paste(substr(i,1,8),'N',sep = ''),'_PSI',sep = '')
  colnames(result_A5SS[[i]])[22] <- paste(paste(substr(i,1,8),'T',sep = ''),'_PSI',sep = '')
  result_A5SS[[i]]$Identity <- NA
  for (number_row in 1:nrow(result_A5SS[[i]])){
    result_A5SS[[i]][number_row,25] <- paste(as.character(result_A5SS[[i]][number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
  }
}
#取出文件的psi的两列以及Identity列，进行循环合并
for(i in names(result_A5SS)){
  tmp_only_psi_dataframe <- result_A5SS[[i]][,c(25,21,22)]
  A5SS_psi_matrix <- merge(A5SS_psi_matrix,tmp_only_psi_dataframe,by = "Identity",all.x  = TRUE)
}
write.table(A5SS_psi_matrix, "C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\A5SS_psi_matrix.txt", sep="\t",quote = FALSE, row.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#读取MXE合并后的文件，添加identity列
MXE <- read.table(file = "C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\head_all_MXE.MATS.JC.txt",header = TRUE)
MXE$Identity <- NA
for (number_row in 1:nrow(MXE)){
  MXE[number_row,27] <- paste(as.character(MXE[number_row,c(3,4,5,6,7,8,9,10,11,12,13)]),collapse = '_')
}
#unique事件
Identity <- unique(sort(MXE$Identity))
MXE_psi_matrix <- data.frame(Identity = Identity)
#循环读取文件
filenames <- c("H9DW0030","H9DW0038","H9DW0039","H9DW0047","H9DW0049","H9DW0052","H9DW0056","H9DW0057","H9DW0064","H9DW0068","H9DW0073","H9DW0074","H9DW0076","H9DW0085","H9DW0092","H9DW0109","H9DW0110","H9DW0116","H9DW0117","H9DW0123","H9DW0125","H9DW0126","H9DW0154","H9DW0169","H9DW0192","H9DW0193","H9DW0194","H9DW0201","H9DW0212","H9DW0219","H9DW0227","H9DW0230","H9DW0250","H9DW0251","H9DW0261","H9DW0262","H9DW0273","H9DW0274","H9DW0280","H9DW0300","H9DW0301","H9DW0312","H9DW0314","H9DW0330","H9DW0338","H9DW0344","H9DW0348","H9DW0351","H9DW0362","H9DW0368","H9DW0377","H9DW0378","H9DW0380","H9DW0383","H9DW0385","H9DW0388","H9DW0079","H9DW0327")
dirpath <- 'C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\MXE'
result_MXE <- list()
for (filename in filenames){
  df <- read.table(paste(dirpath,paste(filename,'_MXE.MATS.JC.txt',sep = ''),sep = "\\"),header = TRUE)
  result_MXE[[paste(filename,'_MXE.MATS.JC.txt',sep = '')]] <- df
}
#将每个文件的inclevel改名为psi并且添加identity列
for(i in names(result_MXE)){
  colnames(result_MXE[[i]])[23] <- paste(paste(substr(i,1,8),'N',sep = ''),'_PSI',sep = '')
  colnames(result_MXE[[i]])[24] <- paste(paste(substr(i,1,8),'T',sep = ''),'_PSI',sep = '')
  result_MXE[[i]]$Identity <- NA
  for (number_row in 1:nrow(result_MXE[[i]])){
    result_MXE[[i]][number_row,27] <- paste(as.character(result_MXE[[i]][number_row,c(3,4,5,6,7,8,9,10,11,12,13)]),collapse = '_')
  }
}
#取出文件的psi的两列以及Identity列，进行循环合并
for(i in names(result_MXE)){
  tmp_only_psi_dataframe <- result_MXE[[i]][,c(27,23,24)]
  MXE_psi_matrix <- merge(MXE_psi_matrix,tmp_only_psi_dataframe,by = "Identity",all.x  = TRUE)
}
write.table(MXE_psi_matrix, "C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\MXE_psi_matrix.txt", sep="\t",quote = FALSE, row.names = FALSE)
#----------------------------------------------------------------------------------------------------------------------------------------------------
#读取RI合并后的文件，添加identity列
RI <- read.table(file = "C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\head_all_RI.MATS.JC.txt",header = TRUE)
RI$Identity <- NA
for (number_row in 1:nrow(RI)){
  RI[number_row,25] <- paste(as.character(RI[number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
}
#unique事件
Identity <- unique(sort(RI$Identity))
RI_psi_matrix <- data.frame(Identity = Identity)
#循环读取文件
filenames <- c("H9DW0030","H9DW0038","H9DW0039","H9DW0047","H9DW0049","H9DW0052","H9DW0056","H9DW0057","H9DW0064","H9DW0068","H9DW0073","H9DW0074","H9DW0076","H9DW0085","H9DW0092","H9DW0109","H9DW0110","H9DW0116","H9DW0117","H9DW0123","H9DW0125","H9DW0126","H9DW0154","H9DW0169","H9DW0192","H9DW0193","H9DW0194","H9DW0201","H9DW0212","H9DW0219","H9DW0227","H9DW0230","H9DW0250","H9DW0251","H9DW0261","H9DW0262","H9DW0273","H9DW0274","H9DW0280","H9DW0300","H9DW0301","H9DW0312","H9DW0314","H9DW0330","H9DW0338","H9DW0344","H9DW0348","H9DW0351","H9DW0362","H9DW0368","H9DW0377","H9DW0378","H9DW0380","H9DW0383","H9DW0385","H9DW0388","H9DW0079","H9DW0327")
dirpath <- 'C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\RI'
result_RI <- list()
for (filename in filenames){
  df <- read.table(paste(dirpath,paste(filename,'_RI.MATS.JC.txt',sep = ''),sep = "\\"),header = TRUE)
  result_RI[[paste(filename,'_RI.MATS.JC.txt',sep = '')]] <- df
}
#将每个文件的inclevel改名为psi并且添加identity列
for(i in names(result_RI)){
  colnames(result_RI[[i]])[21] <- paste(paste(substr(i,1,8),'N',sep = ''),'_PSI',sep = '')
  colnames(result_RI[[i]])[22] <- paste(paste(substr(i,1,8),'T',sep = ''),'_PSI',sep = '')
  result_RI[[i]]$Identity <- NA
  for (number_row in 1:nrow(result_RI[[i]])){
    result_RI[[i]][number_row,25] <- paste(as.character(result_RI[[i]][number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
  }
}
#取出文件的psi的两列以及Identity列，进行循环合并
for(i in names(result_RI)){
  tmp_only_psi_dataframe <- result_RI[[i]][,c(25,21,22)]
  RI_psi_matrix <- merge(RI_psi_matrix,tmp_only_psi_dataframe,by = "Identity",all.x  = TRUE)
}
write.table(RI_psi_matrix, "C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\RI_psi_matrix.txt", sep="\t",quote = FALSE, row.names = FALSE)
#-------------------------------------------------------------------------------------------------------------------------------------------
#读取SE合并后的文件，添加identity列
SE <- read.table(file = "C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\head_all_SE.MATS.JC.txt",header = TRUE)
SE$Identity <- NA
for (number_row in 1:nrow(SE)){
  SE[number_row,25] <- paste(as.character(SE[number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
}
#unique事件
Identity <- unique(sort(SE$Identity))
SE_psi_matrix <- data.frame(Identity = Identity)
#循环读取文件
filenames <- c("H9DW0030","H9DW0038","H9DW0039","H9DW0047","H9DW0049","H9DW0052","H9DW0056","H9DW0057","H9DW0064","H9DW0068","H9DW0073","H9DW0074","H9DW0076","H9DW0085","H9DW0092","H9DW0109","H9DW0110","H9DW0116","H9DW0117","H9DW0123","H9DW0125","H9DW0126","H9DW0154","H9DW0169","H9DW0192","H9DW0193","H9DW0194","H9DW0201","H9DW0212","H9DW0219","H9DW0227","H9DW0230","H9DW0250","H9DW0251","H9DW0261","H9DW0262","H9DW0273","H9DW0274","H9DW0280","H9DW0300","H9DW0301","H9DW0312","H9DW0314","H9DW0330","H9DW0338","H9DW0344","H9DW0348","H9DW0351","H9DW0362","H9DW0368","H9DW0377","H9DW0378","H9DW0380","H9DW0383","H9DW0385","H9DW0388","H9DW0079","H9DW0327")
dirpath <- 'C:\\Users\\xcz\\Desktop\\GCC\\data\\rmats\\SE'
result_SE <- list()
for (filename in filenames){
  df <- read.table(paste(dirpath,paste(filename,'_SE.MATS.JC.txt',sep = ''),sep = "\\"),header = TRUE)
  result_SE[[paste(filename,'_SE.MATS.JC.txt',sep = '')]] <- df
}
#将每个文件的inclevel改名为psi并且添加identity列
for(i in names(result_SE)){
  colnames(result_SE[[i]])[21] <- paste(paste(substr(i,1,8),'N',sep = ''),'_PSI',sep = '')
  colnames(result_SE[[i]])[22] <- paste(paste(substr(i,1,8),'T',sep = ''),'_PSI',sep = '')
  result_SE[[i]]$Identity <- NA
  for (number_row in 1:nrow(result_SE[[i]])){
    result_SE[[i]][number_row,25] <- paste(as.character(result_SE[[i]][number_row,c(3,4,5,6,7,8,9,10,11)]),collapse = '_')
  }
}
#取出文件的psi的两列以及Identity列，进行循环合并
for(i in names(result_SE)){
  tmp_only_psi_dataframe <- result_SE[[i]][,c(25,21,22)]
  SE_psi_matrix <- merge(SE_psi_matrix,tmp_only_psi_dataframe,by = "Identity",all.x  = TRUE)
}
write.table(SE_psi_matrix, "C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\SE_psi_matrix.txt", sep="\t",quote = FALSE, row.names = FALSE)
}
rm(list = ls())
A3SS_psi_matrix <- read.table("C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\A3SS_psi_matrix.txt",header = TRUE)
A5SS_psi_matrix <- read.table("C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\A5SS_psi_matrix.txt",header = TRUE)
MXE_psi_matrix <- read.table("C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\MXE_psi_matrix.txt",header = TRUE)
RI_psi_matrix <- read.table("C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\RI_psi_matrix.txt",header = TRUE)
SE_psi_matrix <- read.table("C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\SE_psi_matrix.txt",header = TRUE)

A3SS_psi_matrix$Filter <- NA
A5SS_psi_matrix$Filter <- NA
MXE_psi_matrix$Filter <- NA
RI_psi_matrix$Filter <- NA
SE_psi_matrix$Filter <- NA

A3SS_psi_matrix_missingcounts <- rowSums(is.na(A3SS_psi_matrix))
A5SS_psi_matrix_missingcounts <- rowSums(is.na(A5SS_psi_matrix))
MXE_psi_matrix_missingcounts <- rowSums(is.na(MXE_psi_matrix))
RI_psi_matrix_missingcounts <- rowSums(is.na(RI_psi_matrix))
SE_psi_matrix_missingcounts <- rowSums(is.na(SE_psi_matrix))

add_label <- function(matrix_name,matrix_missingcounts){
  for(row_number in 1:nrow(matrix_name)){
    if(matrix_missingcounts[row_number] / (ncol(matrix_name) - 2) < 0.2){
      matrix_name[row_number,'Filter'] <- 'pass'
    }else{
      matrix_name[row_number,'Filter'] <- 'fail'
    }
  }
  return(matrix_name)
}

A3SS_psi_matrix <- add_label(A3SS_psi_matrix,A3SS_psi_matrix_missingcounts)
A5SS_psi_matrix <- add_label(A5SS_psi_matrix,A5SS_psi_matrix_missingcounts)
MXE_psi_matrix <- add_label(MXE_psi_matrix,MXE_psi_matrix_missingcounts)
RI_psi_matrix <- add_label(RI_psi_matrix,RI_psi_matrix_missingcounts)
SE_psi_matrix <- add_label(SE_psi_matrix,SE_psi_matrix_missingcounts)

write.table(A3SS_psi_matrix,"C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\label_A3SS_psi_matrix.txt",quote = FALSE,row.names = FALSE,sep = '\t')
write.table(A5SS_psi_matrix,"C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\label_A5SS_psi_matrix.txt",quote = FALSE,row.names = FALSE,sep = '\t')
write.table(MXE_psi_matrix,"C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\label_MXE_psi_matrix.txt",quote = FALSE,row.names = FALSE,sep = '\t')
write.table(RI_psi_matrix,"C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\label_RI_psi_matrix.txt",quote = FALSE,row.names = FALSE,sep = '\t')
write.table(SE_psi_matrix,"C:\\Users\\xcz\\Desktop\\GCC\\psi_matrix\\label_SE_psi_matrix.txt",quote = FALSE,row.names = FALSE,sep = '\t')


A3SS_psi_matrix <- subset(A3SS_psi_matrix,Filter == 'pass')
A5SS_psi_matrix <- subset(A5SS_psi_matrix,Filter == 'pass')
MXE_psi_matrix <- subset(MXE_psi_matrix,Filter == 'pass')
RI_psi_matrix <- subset(RI_psi_matrix,Filter == 'pass')
SE_psi_matrix <- subset(SE_psi_matrix,Filter == 'pass')
