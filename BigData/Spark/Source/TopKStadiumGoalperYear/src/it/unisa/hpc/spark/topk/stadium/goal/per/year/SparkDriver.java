package it.unisa.hpc.spark.topk.stadium.goal.per.year;


import java.util.*;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.*;
import scala.Tuple2;


public class SparkDriver {
        
    public static void main(String[] args) {
        
        String inputFile = args[0];
        String outputPath = args[1];

// Create a configuration object and set the name of the application 
        SparkConf conf = new SparkConf().setAppName("TopKStadiumGoalperYear");

// Create a Spark Context object
        JavaSparkContext sc = new JavaSparkContext(conf);

        // Build an RDD of Strings from the input textual file 
// Each element of the RDD is a line of the input file 
        JavaRDD<String> lines = sc.textFile(inputFile);
        
        String header = lines.first();
        JavaRDD<String> newData = lines.filter(row->!row.equals(header));
        
        JavaPairRDD<String,Integer> stadiumGoalPerMatch = newData.flatMapToPair(row->{
            String[] fields = row.split(",");
            List<Tuple2<String, Integer>> result = new ArrayList<>();
            if(fields[0].equals(args[2]))
                result.add(new Tuple2(fields[3], Integer.valueOf(fields[6])+Integer.valueOf(fields[7])));
            return result.iterator();
        });
        
        JavaPairRDD<String,Integer> totalGoalPerStadium = stadiumGoalPerMatch.reduceByKey((x,y)->x+y);
        List<Tuple2<String, Integer>> topStadiums = totalGoalPerStadium.top(Integer.valueOf(args[3]), new TupleComparator());
        
        JavaPairRDD<String,Integer> topStadiumRDD = sc.parallelizePairs(topStadiums);

        topStadiumRDD.saveAsTextFile(outputPath);
        
        
// Close the Spark Context object
        sc.close();
        
    }
    
}
