options(stringsAsFactors = FALSE)
script_path <- function() {
  z <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  if (length(z)) return(normalizePath(sub('^--file=','',z[1L]),winslash='/',mustWork=TRUE))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile,winslash='/',mustWork=TRUE))
  stop('Could not determine runner path.')
}
test_dir <- dirname(script_path())
root <- normalizePath(file.path(test_dir,'..'),winslash='/',mustWork=TRUE)
results_dir <- file.path(root,'results'); dir.create(results_dir,showWarnings=FALSE,recursive=TRUE)
logs_dir <- file.path(results_dir,'native_logs'); dir.create(logs_dir,showWarnings=FALSE,recursive=TRUE)

rscript <- file.path(R.home('bin'), if (.Platform$OS.type=='windows') 'Rscript.exe' else 'Rscript')
if (!file.exists(rscript)) stop('Could not locate Rscript under R.home().')

tests <- c(
  'TEST_01_MULTITARGET_EXACT.R',
  'TEST_02_ALL_ARCHITECTURE_STATE_CONTRACTS.R',
  'TEST_03_NATIVE_ALL_ARCHITECTURES.R',
  'TEST_04_VERTEBRATE_POINT_PSOCK_OUTPUT_PATCH.R'
)
records <- vector('list',length(tests))
cat('INApest PoF biocontrol/multi-target native validation\n')
cat('R:',R.version.string,'\n')
cat('Platform:',R.version$platform,'\n\n')
for (i in seq_along(tests)) {
  f <- file.path(test_dir,tests[i]); logf <- file.path(logs_dir,sub('\\.R$','.log',tests[i]))
  cat(sprintf('[%d/%d] %s\n',i,length(tests),tests[i]))
  t0 <- proc.time()[['elapsed']]
  status <- system2(rscript,c('--vanilla',shQuote(f)),stdout=logf,stderr=logf)
  elapsed <- proc.time()[['elapsed']] - t0
  ok <- identical(as.integer(status),0L)
  records[[i]] <- data.frame(test=tests[i],status=if(ok)'PASS' else 'FAIL',elapsed_seconds=elapsed,log=basename(logf),stringsAsFactors=FALSE)
  cat(if(ok)'  PASS' else '  FAIL',sprintf('(%.2f s)\n',elapsed))
  if (!ok) {
    cat('  Log:',logf,'\n')
    break
  }
}
rec <- do.call(rbind,records[seq_len(sum(vapply(records,function(x)!is.null(x),logical(1))))])
write.csv(rec,file.path(results_dir,'NATIVE_VALIDATION_BLOCKS.csv'),row.names=FALSE)
overall <- nrow(rec)==length(tests) && all(rec$status=='PASS')
summary <- c(
  'INApest PoF biocontrol/multi-target native validation',
  paste0('Date: ',Sys.Date()),
  paste0('R: ',R.version.string),
  paste0('Platform: ',R.version$platform),
  paste0('Blocks passed: ',sum(rec$status=='PASS'),'/',length(tests)),
  paste0('OVERALL: ',if(overall)'PASS' else 'FAIL')
)
writeLines(summary,file.path(results_dir,'NATIVE_VALIDATION_STATUS.txt'))
cat('\n',paste(summary,collapse='\n'),'\n',sep='')
if (!overall) quit(status=1L)
