package it.unisa.hpc.hadoop.referee.analysis;

import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.io.Text;
import java.io.IOException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.io.LongWritable;

public class RefereeMapper extends Mapper<LongWritable, Text, Text, Text>{
    private String targetNationality; //Nationality taken from input
    private Text refereeName = new Text(); //The name of the Referee/Assistant1/Assistant2
    private Text matchRole = new Text(); //The role of the refereeName as Referee/Assistant1/Assistant2 in a match

    @Override
    protected void setup(Context context) throws IOException, InterruptedException {
        // Retrieve the target nationality from the configuration
        Configuration conf = context.getConfiguration();
        targetNationality = conf.get("nationality");
    }

    @Override
    protected void map(LongWritable key, Text value, Context context) throws IOException, InterruptedException {
        // Convert the input line to a string
        String line = value.toString();

        // Ignore the header
        if (key.get() == 0) {
            return;
        }
        
       // Creation of an array containing the names of Referee, Assistant1 and Assistant2
        String[] fields = line.split(",");  // Split the line into fields using a comma as the delimiter
        String[] referees = {fields[13],fields[14],fields[15]}; // These fields contain info about Referee-Assistant1-Assistant2
        String[] refereeNameFields; //String that will contain all the infos about a single referee
        String refereeNationality; //Variable used for obtaining the Referee/Assistant1/Assistant2 nationality
        
        // Iterate through the referees (Referee, Assistant1, Assistant2)
        for(String referee: referees){
            refereeNameFields = referee.split(" "); //In the single cell of the dataset, name, surname and nationality
                                                                        //of a Referee/Assistant1/Assistant2 are separate with " ",
                                                                        //so we take all the fields one by one
            refereeNationality = refereeNameFields[refereeNameFields.length - 1]; //From the previous array we take the
                                                                                                                           //nationality considering the last element
                                                                                                                           //of the refereeNameFields array
            
            // Check if the referee's nationality matches the target nationality
            if (refereeNationality.contains(targetNationality)) { //We check if the refereeNationality variable taken is equal to
                                                                                           //the targetNationality chosen from input
                refereeName.set(referee); 

                String matchId = fields[17]; //We take the matchId from the right field to associate it with the 
                                                          //referee/assistant1/assistant that ruled that match
                String role = ""; //Initialization of a String role variable to contain the role of the Referee/Assistant1/Assistant2
                
                //Section to check the role
                if (referee.equals(fields[13])) { //If the referee is equal to `fields[13]`, then their role is that of the referee.
                    role = ": Referee";
                } else if (referee.equals(fields[14])) { //If the referee is equal to `fields[14]`, then their role is that of the Assistant1.
                    role = ": Assistant1";
                } else if (referee.equals(fields[15])) { //If the referee is equal to `fields[15]`, then their role is that of the Assistant2..
                    role = ": Assistant2";
                }

                matchRole.set(matchId + role); //Setting of the role of the match
                context.write(refereeName, matchRole); //Writing the role and the referee name into the context
            }
        }
    }
}
