path = "/data/tobiasoh/"
source(sprintf("%sSimulation_Study.R", path))

num_runs = 20
n = 200
load("/data/tobiasoh/truePara.RData")

for (i in 1:num_runs) {
  dataset = make_dataset(n, truePara)
  save(dataset, file=sprintf("/data/tobiasoh/SimStudy/datasets/dataset%d.RData", i))
}