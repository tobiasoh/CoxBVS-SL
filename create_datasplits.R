library(dplyr)

load("./RealData/march8/dataset_basal.RData") #load X, meta

num_runs = 20
training_percentage = 0.9
for (i in 1:num_runs) {
   events_idx = which(meta$vital_status == "Dead")
   censored_idx = dplyr::setdiff(1:nrow(X), events_idx)
   
   train_events = sample(events_idx, size=training_percentage*length(events_idx))
   train_censored = sample(censored_idx, size=training_percentage*length(censored_idx) )
   
   train_indices = sort(c(train_events, train_censored))
   
   #test_indices = append( test_indices, list( dplyr::setdiff(1:nrow(X), train_indices) ) )
   
   test_indices = dplyr::setdiff(1:nrow(X), train_indices)
   save(test_indices, file=sprintf("./RealData/datasplits/test_indices%d.RData", i) ) 
}
