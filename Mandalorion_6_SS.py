import numpy
import sys

content_file=sys.argv[1]
out_path=sys.argv[2]
cutoff=float(sys.argv[3])
genome_file=sys.argv[4]
refine=sys.argv[5]

minimum_read_count=3
minimum_read_coverage=2
splice_site_width=5




def scan_for_best_bin(entry,distance_range,iterator_shift,density_dict,extra_base_cutoff_bottom,extra_base_cutoff_top,peak_areas,chromosome,side):
    best_extra_list=[]
    peak_center=0
    bases=[]
    coverage_area=[]
    best_direction_l=[] 
    best_direction_r=[]


    for x in distance_range:
        extra_list_bases=[]
        extra_list_expression=[]
        extra_list_starts=[]
        extra_list_ends=[]
        direction_l={}
        direction_r={} 
        coverage_set=[]
        bla=0
        for y in distance_range:
            try:
               bla=peak_areas[chromosome][side][entry+x+y]  
            except:
               pass

        if bla==0:  
             highest_y=0
             highest_y_pos=0
             for y in distance_range:
                try: 
 

                    for item in density_dict[entry+x+y]:
                        
                        extra_list_bases.append(item[0])
                        extra_list_expression.append(1)
                        extra_list_starts.append(item[1])
                        extra_list_ends.append(item[2])
                        try:
                            direction_l[item[4]]+=1
                        except:
                            direction_l[item[4]]=1 
                        try:
                            direction_r[item[5]]+=1
                        except:
                            direction_r[item[5]]=1
                        for covered_position in item[3]:

                            coverage_set.append(covered_position)
                except:
                        pass

        if extra_base_cutoff_bottom<=numpy.median(extra_list_bases)<=extra_base_cutoff_top:
            if sum(extra_list_expression)>sum(best_extra_list):
                best_extra_list=extra_list_expression
                peak_center=entry+x

                bases=extra_list_bases
                coverage_area=coverage_set
                best_direction_l=direction_l 
                best_direction_r=direction_r


    return best_extra_list,peak_center,bases,coverage_area,best_direction_l,best_direction_r

def determine_coverage(coverage_area,chromosome,reverse,peak_center,histo_coverage):
    coverage=[]
    coverage.append(0)
    forward=0
    reverse1=0
    forward_missed=0
    reverse_missed=0
    coverage_area2=[]
    coverage_area3=[]
    for covered_position in set(coverage_area):
         coverage_area3.append(covered_position)
         if coverage_area.count(covered_position)>1:
               coverage_area2.append(covered_position)

 
    coverage_area=sorted(coverage_area2,reverse=reverse)
    
    counter=0
    for base_f in coverage_area:
         count=0 
         if reverse==False:
             if base_f>peak_center:
                    count=1
         elif reverse==True:
             if base_f<peak_center:
                    count=1

         if count==1:
             if counter<=3:
                  counter+=1
                  base_f=myround(base_f)
                  try:
                      bla=histo_coverage[chromosome][base_f]

                      coverage.append(histo_coverage[chromosome][base_f])
                  except:
                      pass
             else:
                   break
    coverage=max(coverage)
    return coverage, coverage_area


def read_seq_file(seq_file):

    read_seq={}
        
    length=0
    for line2 in open(seq_file):
        length+=1
    seq_file_open=open(seq_file,'r')
    counter=0
    while counter<length:
        fasta_name=seq_file_open.readline().strip()
        fasta_seq=seq_file_open.readline().strip()
        fasta_name=fasta_name[1:] 
        read_seq[fasta_name]=fasta_seq
        counter+=2
    return read_seq


def reverse_complement(sequence):
  Seq=''
  complement = {'a':'T','c':'G','g':'C','t':'A','n':'N','A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq



def myround(x, base=10):
    return int(base * round(float(x)/base))


def find_peaks(density_dict,out,Peaks,reverse,cutoff,extra_base_cutoff_bottom,extra_base_cutoff_top,histo_coverage,side,peak_areas,chromosome):

    if reverse==False:
        distance_range=range(-splice_site_width,splice_site_width,1)
        iterator_shift=1
    if reverse==True:
        distance_range=range(splice_site_width,-splice_site_width,-1)
        iterator_shift=-1


    entry_list=[]
    for entry in density_dict:
      entry_list.append([entry,density_dict[entry]])

    for entry,density in sorted(entry_list,key=lambda x: sum(numpy.array(x[1])[:,2]),reverse=True):
        if len(density)>=minimum_read_count:

      
            try:
                bla=peak_areas[chromosome][side][entry]
            except:

                best_extra_list,peak_center,bases,coverage_area,best_direction_l,best_direction_r=scan_for_best_bin(entry,distance_range,iterator_shift,density_dict,extra_base_cutoff_bottom,extra_base_cutoff_top,peak_areas,chromosome,side)

                coverage,coverage_area=determine_coverage(coverage_area,chromosome,reverse,peak_center,histo_coverage)


                if coverage>0:
                    print(chromosome+'\t'+str(peak_center-1)+'\t'+str(peak_center+1)+'\t'+str(round(sum(best_extra_list)/coverage,3))) 
                    if round(sum(best_extra_list)/coverage,3)>cutoff:
                         try:
                             Left_TSS=best_direction_l['TSS']
                         except:
                             Left_TSS=0
                         try:
                             Left_TES=best_direction_l['TES']
                         except:
                             Left_TES=0
                         try:
                             Right_TSS=best_direction_r['TSS']
                         except:
                             Right_TSS=0
                         try:
                             Right_TES=best_direction_r['TES']
                         except:
                             Right_TES=0
                         Left_to_Right=Left_TSS+Right_TES
                         Right_to_Left=Left_TES+Right_TSS
                         Type='-'
          
                         if Left_to_Right>Right_to_Left:
                             if reverse==True:
                                 Type='5'
                             elif reverse==False:
                                 Type='3'
                         if Left_to_Right<Right_to_Left:
                             if reverse==True:
                                 Type='3'
                             elif reverse==False:
                                 Type='5'

                         if Type!='-':
                             Peaks+=1

                             out.write(chromosome+'\t'+str(peak_center-splice_site_width)+'\t'+str(peak_center+splice_site_width)+'\t'+str(Type)+side+str(Peaks)+'_'+str(peak_center-splice_site_width)+'_'+str(peak_center+splice_site_width)+'_'+str(round(sum(best_extra_list)/coverage,3))+'\t'+str(Peaks)+'\n')
                             


                             for base in range(peak_center-splice_site_width,peak_center+splice_site_width,1):                             
                                 peak_areas[chromosome][side][base]=1

        else:
            break                    

    return Peaks,peak_areas





def collect_reads(content_file):
    histo_left_bases={}
    histo_right_bases={}
    chromosome_list_left=set()
    chromosome_list_right=set()
    histo_coverage={}

    extra_base_cutoff_top=5
    extra_base_cutoff_bottom=0
    for line in open(content_file):
#        print(line.strip())

        total=0
        b=line.strip().split('\t')
        infile=b[0]

        length=0

        for line in open(infile):
          total+=1
          a=line.strip().split('\t')
          chromosome=a[13]
          
          if not histo_coverage.get(chromosome):
              histo_coverage[chromosome]={}


          score=int(a[0])
          direction=a[8]
          name=a[9]
          matches=name.split('_')
          coverage=int(matches[3])
          if coverage>=minimum_read_coverage:
            matches=name.split('_')
            length=int(matches[5].split('|')[0])
            begin=int(a[15])
            span=int(a[16])
            blocksizes=a[18].split(',')[:-1]
            blockstarts=a[20].split(',')[:-1]
            readstarts=a[19].split(',')[:-1]

            if direction=='+':
                start_seq='S'
                end_seq='E'
                left_match='TSS'
                right_match='TES'

            else:
                start_seq='E'
                end_seq='S'
                left_match='TES'
                right_match='TSS'



            coverage_set=set()
            previous_blocksize=-1
            previous_start=-1
            previous_blockend=100000000000000
            intron=0
            indel=0
            low_bounds=[]
            up_bounds=[]
            indel1=0

            aligned_bases=0

            for x in range(0,len(blocksizes),1):
                blockstart=int(blockstarts[x])
                blocksize=int(blocksizes[x])
                readstart=int(readstarts[x])
                aligned_bases+=blocksize
                blockend=blockstart+blocksize
                if blocksize>10:


                   for y in range(0,blocksize,10):

                       rounded=myround(blockstart+y)
                       coverage_set.add(rounded)

                   for yy in range(y,blocksize,1):
                       rounded=myround(blockstart+yy)
                       coverage_set.add(rounded)
    





                   if previous_start==-1:
                       previous_start=blockstart
                       min_length=10
                   else:
                       min_length=10
 
                   if blockstart-previous_blockend>20:
                       previous_start=blockstart
 
                   if blockend-previous_start>min_length:
                       if intron==1:
                          up_bounds.append([previous_start,indel1,blockend])
                          low_bounds.append([remember_blockend,remember_indel1,remember_start])
                          intron=0
                       else:
                          try:
                             next_blockstart=int(blockstarts[x+1])
                             next_blocksize=int(blocksizes[x+1])
                             next_readstart=int(readstarts[x+1])

                             insert=next_blockstart-blockend
                             if insert>50:
                                 indel1=next_readstart-(readstart+blocksize) 
                                 remember_blockend=blockend
                                 remember_indel1=indel1
                                 remember_start=previous_start
                                 intron=1
                                 previous_start=next_blockstart 
                                 blockend=next_blockstart  
                          except:
                              pass

                   previous_blockend=blockend 

            for rounded in coverage_set:           
                try:
                    histo_coverage[chromosome][rounded]+=1
                except:
                    histo_coverage[chromosome][rounded]=1

            ratio=aligned_bases/length

            if ratio>0.70:
  
                for low_bound,indel1,blockend in low_bounds:  

                    chromosome_list_left.add(chromosome)
                    if not histo_left_bases.get(chromosome):
                        histo_left_bases[chromosome]={}  
                    if not histo_left_bases[chromosome].get(low_bound):
                        histo_left_bases[chromosome][low_bound]=[]
                    histo_left_bases[chromosome][low_bound].append([indel1,begin,span,coverage_set,left_match,right_match])



                for up_bound,indel1,blockend in up_bounds:   

                    chromosome_list_right.add(chromosome)            
                    if not histo_right_bases.get(chromosome):
                        histo_right_bases[chromosome]={}  
                    if not histo_right_bases[chromosome].get(up_bound):
                        histo_right_bases[chromosome][up_bound]=[]
                    histo_right_bases[chromosome][up_bound].append([indel1,begin,span,coverage_set,left_match,right_match])


    chromosome_list=chromosome_list_left & chromosome_list_right
    return histo_left_bases, histo_right_bases,chromosome_list,histo_coverage


def parse_genome(input_file,left_bounds,right_bounds):


    gene_dict={} 
    for line in open(input_file):
        a=line.strip().split('\t')
        if len(a)>7:
             if a[2]=='exon':
                 if not gene_dict.get(a[8].split('; transcript_id "')[1].split('"')[0]):
                     gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]]=[]

                 gene_dict[a[8].split('; transcript_id "')[1].split('"')[0]].append((a[0],a[3],a[4],a[6]))

    read_list=[]
    for transcript_id in gene_dict:
        transcript_data=gene_dict[transcript_id]

        chromosome=transcript_data[0][0]
        try:
            bla=right_bounds[chromosome]
        except:
            left_bounds[chromosome]={}
            right_bounds[chromosome]={}
            left_bounds[chromosome]['5']=[]
            right_bounds[chromosome]['5']=[]
            left_bounds[chromosome]['3']=[]
            right_bounds[chromosome]['3']=[]

        start=sorted(transcript_data,key=lambda x: int(x[1]))[0][1]
        end=sorted(transcript_data,key=lambda x: int(x[2]),reverse=True)[0][2]

        for entry in transcript_data:
              if entry[1]!=start:
                  if entry[3]=='+':
                      right_bounds[chromosome]['3'].append(int(entry[1])-1)  
                  elif entry[3]=='-':
                      right_bounds[chromosome]['5'].append(int(entry[1])-1)  
              if entry[2]!=end:
                  if entry[3]=='+':
                      left_bounds[chromosome]['5'].append(int(entry[2]))  
                  if entry[3]=='-':
                      left_bounds[chromosome]['3'].append(int(entry[2]))  
     
    return left_bounds,right_bounds

def make_genome_bins(bounds,side,Peaks,chromosome,peak_areas):
    for type1 in ['5','3']:
        covered={}
        position_list=sorted(bounds[type1],key=int)
        for index1 in range(0,len(position_list),1):

            try:
                bla=covered[index1]
            except:
                sub_list=[]
                sub_list.append(position_list[index1])
                for index2 in range(index1,len(position_list),1):
                    if position_list[index2]-max(sub_list)<=splice_site_width:
                        sub_list.append(position_list[index2])
                        covered[index2]=1
                    else:
                        break
                single=0                   
                if len(sub_list)>1:
                    splice_distances=[]
                    for splice_pos in range(0,len(sub_list)-1,1):
                        splice_distances.append(int(sub_list[splice_pos+1])-int(sub_list[splice_pos]))
                    if min(splice_distances)>3:
 
                        for x in range(0,len(sub_list),1):
                            if x!=0:
                                start=int(sub_list[x]-((sub_list[x]-sub_list[x-1])/2))
                            else:
                                start=int(sub_list[x])-1
                            if x!=len(sub_list)-1:
                                end=int(sub_list[x]+((sub_list[x+1]-sub_list[x])/2)) 
                            else:
                                end=int(sub_list[x])+1

                            out.write(chromosome+'\t'+str(start)+'\t'+str(end)+'\t'+type1+side+str(Peaks)+'_'+str(start)+'_'+str(end)+'_A'+'\t'+str(Peaks)+'\n')
                            for base in range(start,end,1):                             
                                peak_areas[chromosome][side][base]=1
   
                            Peaks+=1
                    else:
                         single=1
                else:
                    single=1              
                if single==1:
                    start=min(sub_list)-splice_site_width
                    end=max(sub_list)+splice_site_width
                    out.write(chromosome+'\t'+str(start)+'\t'+str(end)+'\t'+type1+side+str(Peaks)+'_'+str(start)+'_'+str(end)+'_A'+'\t'+str(Peaks)+'\n')
                    for base in range(start-1,end+1,1):                             
                        peak_areas[chromosome][side][base]=1
                    Peaks+=1                    
                    
    return Peaks,peak_areas

left_bounds={}
right_bounds={}

left_bounds,right_bounds=parse_genome(genome_file,left_bounds,right_bounds)

Left_Peaks=0
Right_Peaks=0


histo_left_bases, histo_right_bases,chromosome_list,histo_coverage=collect_reads(content_file)
out=open(out_path+'/SS.bed','w')


peak_areas={}
print(chromosome_list)
for chromosome in chromosome_list:
#for chromosome in ['chr16']:
    peak_areas[chromosome]={}
    peak_areas[chromosome]['l']={}
    peak_areas[chromosome]['r']={}
    print(chromosome)
    if 'g' in refine:

        Left_Peaks_old=Left_Peaks
        Right_Peaks_old=Right_Peaks
        Left_Peaks,peak_areas=make_genome_bins(left_bounds[chromosome],'l',Left_Peaks,chromosome,peak_areas)
        Right_Peaks,peak_areas=make_genome_bins(right_bounds[chromosome],'r',Right_Peaks,chromosome,peak_areas)
        print('Annotation-Based', Left_Peaks-Left_Peaks_old,Right_Peaks-Right_Peaks_old)
    Left_Peaks_old=Left_Peaks
    Right_Peaks_old=Right_Peaks
    Left_Peaks,peak_areas=find_peaks(histo_left_bases[chromosome],out,Left_Peaks,True,cutoff,0,5,histo_coverage,'l',peak_areas,chromosome)
    Right_Peaks,peak_areas=find_peaks(histo_right_bases[chromosome],out,Right_Peaks,False,cutoff,0,5,histo_coverage,'r',peak_areas,chromosome)
    print('Read-Based', Left_Peaks-Left_Peaks_old,Right_Peaks-Right_Peaks_old)


print(Left_Peaks)
print(Right_Peaks)












    
