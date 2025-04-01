package it.unisa.hpc.hadoop.referee.analysis;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;


/*
* Driver to manage the Job "RefereeIndex" used to obtain all the referee names of a certain * nationality. 
* The nationality is taken from the command line.
*/
public class RefereeDriver {

    public static void main(String[] args) throws Exception {
        // Check if the correct number of command line arguments is provided
        if (args.length != 3) {
            System.err.println("Usage: RefereeAnalysisHadoop <input path> <output path> <nationality>");
            System.exit(-1);
        }

        // Set up the Hadoop configuration
        Configuration conf = new Configuration();
        conf.set("nationality", args[2]);  // Set the nationality from user input

        // Create a new MapReduce job named "RefereeIndex"
        Job job = Job.getInstance(conf, "RefereeIndex");

        // Set the classes for the MapReduce job
        job.setJarByClass(RefereeDriver.class);
        job.setMapperClass(RefereeMapper.class);
        job.setReducerClass(RefereeReducer.class);

        // Set the output key and value types
        job.setOutputKeyClass(Text.class);
        job.setOutputValueClass(Text.class);

        // Set the input and output paths
        FileInputFormat.addInputPath(job, new Path(args[0]));
        FileOutputFormat.setOutputPath(job, new Path(args[1]));

        // Execute the job and wait for its completion
        System.exit(job.waitForCompletion(true) ? 0 : 1);
    }
}
