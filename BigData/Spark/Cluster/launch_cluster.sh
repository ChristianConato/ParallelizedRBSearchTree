/opt/bitnami/spark/bin/spark-submit \
  --class it.unisa.hpc.spark.topk.stadium.goal.per.year.SparkDriver \
  --master spark://spark-master:7077 \
  --deploy-mode client \
  --supervise \
  --executor-memory 1G \
  ./TopKStadiumGoalperYear.jar \
  ./input ./outputcluster 1958 9
