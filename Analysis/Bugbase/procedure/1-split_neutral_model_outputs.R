# Print commands to command line
options(echo=TRUE)

# Get user input and assign to variables
args <- commandArgs(trailingOnly=TRUE)

neutral_output <- args[1]
neutral_table_below <- args[2]
neutral_table_neutral <- args[3]
neutral_table_above <- args[4]

neutral_all = read.csv(file=neutral_output, header = T)

neutral_below = neutral_all[neutral_all$model=="below",]
neutral_neutral = neutral_all[neutral_all$model=="neutral",]
neutral_above = neutral_all[neutral_all$model=="above",]

write.csv(neutral_below, file=neutral_table_below, row.names=F)
write.csv(neutral_neutral, file=neutral_table_neutral, row.names=F)
write.csv(neutral_above, file=neutral_table_above, row.names=F)