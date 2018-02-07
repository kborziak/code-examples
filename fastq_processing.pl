#!/usr/bin/perl

$files = {
#  'SRR042289' => {'type' => 'sanger', 'organism' => 'human', 'adapter' => 'illumina'}, 
#  'SRR042290' => {'type' => 'sanger', 'organism' => 'human', 'adapter' => 'illumina'}, 
#  'SRR042291' => {'type' => 'sanger', 'organism' => 'human', 'adapter' => 'illumina'}, 
#  'DRR002473' => {'type' => 'illumina', 'organism' => 'mouse', 'adapter' => 'illumina'}, 
#  'DRR002474' => {'type' => 'illumina', 'organism' => 'mouse', 'adapter' => 'illumina'}, 
#  'dmel_130807' => {'type' => 'sanger', 'organism' => 'fly', 'adapter' => 'none'}, 
#  'GSM280085' => {'type' => 'sanger', 'organism' => 'fly', 'adapter' => 'other'}, 
#  'SRR513183' => {'type' => 'sanger', 'organism' => 'mouse', 'adapter' => 'illumina_v1.5'}, 
#  'SRR513184' => {'type' => 'sanger', 'organism' => 'mouse', 'adapter' => 'illumina_v1.5'}, 
#  'SRR1593702' => {'type' => 'sanger', 'organism' => 'fly', 'adapter' => 'illumina_v1.9'}, 
#  'SRR1593704' => {'type' => 'sanger', 'organism' => 'fly', 'adapter' => 'illumina_v1.9'}, 
#  'SRR1593705' => {'type' => 'sanger', 'organism' => 'fly', 'adapter' => 'illumina_v1.9'}, 
#  'oocyte1' => {'type' => 'sanger', 'organism' => 'fly', 'adapter' => 'illumina_v1.9'}, 
#  'oocyte_2' => {'type' => 'sanger', 'organism' => 'fly', 'adapter' => 'none'}, 
   'SRR1272328' => {'type' => 'illumina', 'organism' => 'fly', 'adapter' => 'illumina_v1.5'}, 
};
$adapters = {'illumina' => 'TCGTATGCCGTCTTCTGCTTG', 'other' => 'CTGTAGGCACCATCAA', 'illumina_v1.5' => 'ATCTCGTATGCCGTCTTCTGCTTG', 'illumina_v1.9' => 'TGGAATTCTCGGGTGCCAAGGAACTCCA'};
$genomes = {
  'human' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/ref_GRCh37.p13', 
  'mouse' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/ref_GRCm38.p1', 
  'fly' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/dmel-all-chromosome-r5.54', 
};
$gffs = {
  'human' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/ref_GRCh37.p13_top_level.gff3',
  'mouse' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/ref_GRCm38.p1_top_level.gff3',
  'fly' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/dmel-all-r5.54.gff', 
};
$mirbase = {
  'human' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/mirbase_hsa.gff3',
  'mouse' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/mirbase_mmu.gff3',
  'fly' => '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/genomes/mirbase_dme.gff3', 
};

$path = '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/sra_data';
#$path = '/media/Ext1/analysis/bowtie';
$out_path = '/home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/data';
#$out_path = '/media/Ext1/analysis/bowtie';
$storage_path = '/media/Ext1/analysis/bowtie';
$processing = 'yes';
$bowtie = 'no';
$analysis = 'none';


foreach $id (keys %$files) {
  print "Running $id\n";
  if ($processing eq 'yes') {
# Trim adapter sequences
    print "Adapter trimming\n";
    if ($files->{$id}->{'adapter'} ne 'none') {
#      system ("cutadapt -f fastq -O 1 -a ".$adapters->{$files->{$id}->{'adapter'}}." -o ".$path."/".$id.".3P_trim ".$path."/".$id.".fastq");
    }
    else {
#      system ("cp ".$path."/".$id.".fastq ".$path."/".$id.".3P_trim");
    }
# Trim low quality sequences
    print "Trimming low quality sequences\n";
    if ($files->{$id}->{'type'} eq 'sanger') {
      system ("btrim64 -q -w 5 -a 20 -l 10 -S -t ".$path."/".$id.".3P_trim -o ".$path."/".$id.".3P_trim.b_trim");
    }
    else {
      system ("btrim64 -q -w 5 -a 20 -l 10 -t ".$path."/".$id.".3P_trim -o ".$path."/".$id.".3P_trim.b_trim");
    }
# Remove short trimmed reads
    print "Removing short reads\n";
    system ("/home/kborziak/Documents/epigenetics/miRNA/sra_quality_trim.pl ".$files->{$id}->{'type'}." 10 10 ".$path."/".$id.".3P_trim.b_trim");
# Rename output file
    system ("mv ".$path."/".$id.".3P_trim.b_trim.q_trim ".$path."/".$id.".trim");
    system ("mv ".$path."/".$id.".3P_trim.b_trim.unique ".$path."/".$id.".unique");
  }
  
#  print $id."\n";
  if ($bowtie eq 'yes') {
    print "Running bowtie\n";
    system ("bowtie -v 1 -a --best --strata -m 10 -p 6 ".$genomes->{$files->{$id}->{'organism'}}." ".$path."/".$id.".unique ".$storage_path."/".$id.".bowtie --max ".$storage_path."/".$id.".bowtie_max --un ".$storage_path."/".$id.".bowtie_un");
  }
  #  Bowtie2 mapping
  #  /usr/local/src/bowtie2-2.2.4/bowtie2 -q --phred33 -N 0 --no-1mm-upfront --end-to-end --very-fast -a -p 1 -x dme_hairpin_mature -S oocyte_1.sam -U /home/kborziak/Documents/epigenetics/miRNA/miRNA-seq/sra_data/oocyte1.trim
  
  if ($analysis ne 'none') {
    print "Running analysis\n";
    if ($analysis eq 'all') {
      system ("/home/kborziak/Documents/epigenetics/miRNA/gff_2_gene_3.pl ".$storage_path."/".$id.".bowtie ".$gffs->{$files->{$id}->{'organism'}}." ".$out_path."/".$id.".genes ".$files->{$id}->{'organism'}." not_all");
      system ("/home/kborziak/Documents/epigenetics/miRNA/gene_totals.pl ".$id." ".$files->{$id}->{'organism'});
    }
    elsif ($analysis eq 'miRNA') {
      system ("/home/kborziak/Documents/epigenetics/miRNA/gff_2_gene_mirbase.pl ".$storage_path."/".$id.".bowtie ".$mirbase->{$files->{$id}->{'organism'}}." ".$out_path."/".$id.".miRNA ".$files->{$id}->{'organism'}." ".$id);
    }
  }
}

