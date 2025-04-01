package it.unisa.hpc.spark.topk.stadium.goal.per.year;

import java.util.Comparator;
import scala.Serializable;
import scala.Tuple2;

public class TupleComparator implements Comparator<Tuple2<String, Integer>>, Serializable  {

    @Override
    public int compare(Tuple2<String, Integer> tuple1, Tuple2<String, Integer> tuple2) {
        Integer result = tuple1._2().compareTo(tuple2._2());
        if (result.equals(0))
            result=tuple2._1().compareTo(tuple1._1());
        return result;
    }   
}
