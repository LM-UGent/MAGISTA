import sys

def main(args):
  is_header=True
  with open(args[1]) as f:
    for line in f:
      if is_header:
        is_header=False
        print("bin,best_match,largest_assigned_completeness,purity,F1,shannon_entropy,Quality,size,weighed_completeness,weighed_purity,weighed_F1,distance_to_match,bin_class")
      else:
        parts = line.strip().split(',')
        binid = parts[0]
        best_match = parts[4]
        largest_assigned_completeness = float(parts[2])/float(parts[3])*100.0
        purity = float(parts[2])/float(parts[1])*100.0
        F1 = 2*largest_assigned_completeness*purity/(largest_assigned_completeness+purity)
        shannon_entropy = -1
        Quality = -1
        size = float(parts[1])
        weighed_completeness = -1
        weighed_purity = -1
        weighed_F1 = -1
        D2M = parts[5]
        binclass = parts[6]
        print(','.join([str(binid),str(best_match),str(largest_assigned_completeness),str(purity),str(F1),str(shannon_entropy),str(Quality),str(size),str(weighed_completeness),str(weighed_purity),str(weighed_F1),str(D2M),str(binclass)]))


if __name__ == "__main__":
  if(len(sys.argv)!=2):
    sys.stderr.write("Usage: python "+sys.argv[0]+" <input>\n")
  else:
    main(sys.argv)

