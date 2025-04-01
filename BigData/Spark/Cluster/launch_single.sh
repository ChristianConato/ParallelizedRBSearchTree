/opt/bitnami/spark/bin/spark-submit \
  --class it.unisa.hpc.spark.topk.stadium.goal.per.year.SparkDriver \
  --master local \
  ./TopKStadiumGoalperYear.jar \
  ./input ./outputsingle 1934 4