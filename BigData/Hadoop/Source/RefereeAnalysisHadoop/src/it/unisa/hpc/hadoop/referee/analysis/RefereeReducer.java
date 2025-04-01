package it.unisa.hpc.hadoop.referee.analysis;

import java.io.IOException;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Reducer;

  public class RefereeReducer extends Reducer<Text, Text, Text, Text> {
        private Text result = new Text();

        @Override
        public void reduce(Text key, Iterable<Text> values, Context context) throws IOException, InterruptedException {
            StringBuilder matches = new StringBuilder("\n");
            
            // Concatena tutte le partite arbitrate dall'arbitro
            for (Text val : values) {
                matches.append(val.toString()).append("\n");
            }

            result.set(matches.toString());
            context.write(key, result);
        }
  }
