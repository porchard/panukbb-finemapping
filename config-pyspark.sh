HAIL_HOME=$(pip3 show hail | grep Location | awk -F' ' '{print $2 "/hail"}')
export PYSPARK_SUBMIT_ARGS=" \
  --jars $HAIL_HOME/backend/hail-all-spark.jar \
  --conf spark.driver.memory=5G \
  --conf spark.ui.enabled=false \
  --conf spark.executor.memory=30G \
  --conf spark.driver.extraClassPath=$HAIL_HOME/hail-all-spark.jar \
  --conf spark.executor.extraClassPath=./hail-all-spark.jar \
  --conf spark.serializer=org.apache.spark.serializer.KryoSerializer \
  --conf spark.kryo.registrator=is.hail.kryo.HailKryoRegistrator pyspark-shell"
