runsiamcat <- function(featTable, metaTable, fileName, case, ml, norm, cutoff.p){
  dim(featTable)
  fileName=paste(fileName, ml, norm, cutoff.p)
  print(fileName)
  
  # create SIAMCAT object and classify
  siamcat <- siamcat(feat=featTable, 
                     meta=metaTable, 
                     label="Event", 
                     case=case)
  # abundance and prevelance filtering
  siamcat <- filter.features(siamcat,
                             filter.method = 'abundance',
                             cutoff=0.0001,
                             verbose=3)
  # 
   siamcat <- filter.features(siamcat,
                              filter.method = 'prevalence',
                              cutoff = cutoff.p,
                              feature.type = 'filtered',
                              verbose=3)

  # normalize with log.clr
  siamcat <- normalize.features(siamcat, 
                                norm.method = norm, 
                                feature.type = 'filtered',
                                norm.param = list(log.n0=1e-05, sd.min.q=0.1))
  
  # compute associations 
  siamcat <- check.associations(siamcat, feature.type = 'normalized',
                                detect.lim = 10^-5, 
                                plot.type = "quantile.box",
                                fn.plot = paste0(PARAM$folder.results,
                                          Sys.Date(), fileName, 'assoc.plot.pdf'))

  # train model
  siamcat <- create.data.split(siamcat, num.folds =10, num.resample = 10)  
  siamcat <- train.model(siamcat, method = ml, verbose = 3)
  siamcat <- make.predictions(siamcat)
  siamcat <- evaluate.predictions(siamcat)    
  print(siamcat@eval_data$auroc)
  # evaluation plot
  model.evaluation.plot(siamcat, fn.plot = paste0(PARAM$folder.results, Sys.Date(), '.',
                                                  fileName, 'eval.plot.pdf'))
  # interpretation plot
  model.interpretation.plot(siamcat, fn.plot = paste0(PARAM$folder.results, 
                                                      Sys.Date(), '.', 
                                                      fileName,
                                                      'interpret.plot.pdf'),
                            consens.thres = 0.5,
                            detect.lim = 1e-05,
                            heatmap.type = 'zscore')
  
  # save siamcat object
  save(siamcat, file = paste0(PARAM$folder.results, Sys.Date(), '.', fileName, 'siamcat.Rdata'))
  return(siamcat) 
}
