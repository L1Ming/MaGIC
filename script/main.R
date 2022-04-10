ss = c(
  'MafProcess.R',
  'Optimal.R',
  'CatcherVs.R',
  'Capturing.R',
  'Grand_plot.R',
  'Pan_unicox.R',
  'NPC_processing.R',
  'NvsT.R'
  )

for(s in ss){
  runScript = function(s){
    rm(list = ls())
    message(rep(c('\n', '-', '=', s, '=', '-'), times = c(1, 10, 5, 1, 5, 10)))
    message('Doing...\r', appendLF = F)
    a = Sys.time()
    source(paste0('script/', s), encoding = 'UTF-8')
    message('Done!       ')
    b = Sys.time()
    message("Time Ellapse: ",b-a)
  }
  
  runScript(s)
}
