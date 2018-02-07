#!/usr/bin/perl

use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use Scalar::Util qw(looks_like_number);

$pwd = '/media/kborziak/e03f17a9-e8c0-46e5-943e-7a1ba4082557/dro_testis/14120R/Fastq';
$pwd2 = '/media/kborziak/e03f17a9-e8c0-46e5-943e-7a1ba4082557/dro_testis/fasta';
$files = {
  '14120X1_170412_1_2' =>  {1 => '14120X1_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X1_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X1_170414_1_2' =>  {1 => '14120X1_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X1_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X2_170412_1_2' =>  {1 => '14120X2_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X2_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X2_170414_1_2' =>  {1 => '14120X2_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X2_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X3_170412_1_2' =>  {1 => '14120X3_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X3_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X3_170414_1_2' =>  {1 => '14120X3_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X3_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X4_170412_1_2' =>  {1 => '14120X4_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X4_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X4_170414_1_2' =>  {1 => '14120X4_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X4_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X5_170412_1_2' =>  {1 => '14120X5_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X5_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X5_170414_1_2' =>  {1 => '14120X5_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X5_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X6_170412_1_2' =>  {1 => '14120X6_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X6_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X6_170414_1_2' =>  {1 => '14120X6_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X6_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X7_170412_1_2' =>  {1 => '14120X7_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X7_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X7_170414_1_2' =>  {1 => '14120X7_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X7_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X8_170412_1_2' =>  {1 => '14120X8_170412_D00294_0311_ACAJ3TANXX_6_1.txt',   2 => '14120X8_170412_D00294_0311_ACAJ3TANXX_6_2.txt',   }, 
  '14120X8_170414_1_2' =>  {1 => '14120X8_170414_D00454R_0145_ACAL3CANXX_2_1.txt',  2 => '14120X8_170414_D00454R_0145_ACAL3CANXX_2_2.txt',  }, 
  '14120X11_170412_1_2' => {1 => '14120X11_170412_D00294_0311_ACAJ3TANXX_6_1.txt',  2 => '14120X11_170412_D00294_0311_ACAJ3TANXX_6_2.txt',  }, 
  '14120X11_170414_1_2' => {1 => '14120X11_170414_D00454R_0145_ACAL3CANXX_2_1.txt', 2 => '14120X11_170414_D00454R_0145_ACAL3CANXX_2_2.txt', }, 
  '14120X12_170412_1_2' => {1 => '14120X12_170412_D00294_0311_ACAJ3TANXX_6_1.txt',  2 => '14120X12_170412_D00294_0311_ACAJ3TANXX_6_2.txt',  }, 
  '14120X12_170414_1_2' => {1 => '14120X12_170414_D00454R_0145_ACAL3CANXX_2_1.txt', 2 => '14120X12_170414_D00454R_0145_ACAL3CANXX_2_2.txt', }, 
};

# Tophap
foreach $out (keys $files) {
  $cmd = "mkdir $pwd/tophat_$out; mv $pwd/".$files->{$out}->{1}." $pwd/tophat_$out/.; mv $pwd/".$files->{$out}->{2}." $pwd/tophat_$out/.; cd $pwd/tophat_$out; /usr/local/src/tophat-2.1.1.Linux_x86_64/tophat -p 10 -G $pwd2/dmel-all-r6.16.gtf $pwd2/dmel-all-r6.16 ".$files->{$out}->{1}." ".$files->{$out}->{2}." 2>&1";
	$output = `$cmd`;
  print $output."\n";
}


# Cufflinks
foreach $out (keys $files) {
  $cmd = "cd $pwd/tophat_$out; /usr/local/src/cufflinks-2.2.1.Linux_x86_64/cufflinks -p 5 -g $pwd2/dmel-all-r6.16.gtf -b $pwd2/dmel-all-r6.16.fa -u $pwd/tophat_$out/tophat_out/accepted_hits.bam";
	$output = `$cmd`;
  print $output."\n";
}

# Cuffmerge
# assembly_GTF_list.txt lists location of transcripts.gtf outputs from cufflinks
$cmd = "cd $pwd; /usr/local/src/cufflinks-2.2.1.Linux_x86_64/cuffmerge -p 5 -g $pwd2/dmel-all-r6.16.gtf -s $pwd2/dmel-all-r6.16.fa assembly_GTF_list.txt";
$output = `$cmd`;
print $output."\n";

# Cuffquant
foreach $out (keys $files) {
  $cmd = "cd $pwd/tophat_$out; /usr/local/src/cufflinks-2.2.1.Linux_x86_64/cuffquant -p 5 -b $pwd2/dmel-all-r6.16.fa -u $pwd/merged_asm/merged.gtf $pwd/tophat_$out/tophat_out/accepted_hits.bam";
	$output = `$cmd`;
  print $output."\n";
}

#Cuffdiff
#$cmd = "/usr/local/src/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 5 merged_asm/merged.gtf tophat_14120X1_170412_1_2/abundances.cxb,tophat_14120X1_170414_1_2/abundances.cxb,tophat_14120X2_170412_1_2/abundances.cxb,tophat_14120X2_170414_1_2/abundances.cxb,tophat_14120X3_170412_1_2/abundances.cxb,tophat_14120X3_170414_1_2/abundances.cxb,tophat_14120X4_170412_1_2/abundances.cxb,tophat_14120X4_170414_1_2/abundances.cxb,tophat_14120X5_170412_1_2/abundances.cxb,tophat_14120X5_170414_1_2/abundances.cxb,tophat_14120X6_170412_1_2/abundances.cxb,tophat_14120X6_170414_1_2/abundances.cxb tophat_14120X7_170412_1_2/abundances.cxb,tophat_14120X7_170414_1_2/abundances.cxb,tophat_14120X8_170412_1_2/abundances.cxb,tophat_14120X8_170414_1_2/abundances.cxb,tophat_14120X11_170412_1_2/abundances.cxb,tophat_14120X11_170414_1_2/abundances.cxb,tophat_14120X12_170412_1_2/abundances.cxb,tophat_14120X12_170414_1_2/abundances.cxb -o cuffdiff/";
$cmd = "/usr/local/src/cufflinks-2.2.1.Linux_x86_64/cuffdiff -p 5 merged_asm/merged.gtf tophat_14120X2_170412_1_2/abundances.cxb,tophat_14120X2_170414_1_2/abundances.cxb,tophat_14120X3_170412_1_2/abundances.cxb,tophat_14120X3_170414_1_2/abundances.cxb tophat_14120X4_170412_1_2/abundances.cxb,tophat_14120X4_170414_1_2/abundances.cxb,tophat_14120X5_170412_1_2/abundances.cxb,tophat_14120X5_170414_1_2/abundances.cxb,tophat_14120X6_170412_1_2/abundances.cxb,tophat_14120X6_170414_1_2/abundances.cxb tophat_14120X7_170412_1_2/abundances.cxb,tophat_14120X7_170414_1_2/abundances.cxb,tophat_14120X8_170412_1_2/abundances.cxb,tophat_14120X8_170414_1_2/abundances.cxb tophat_14120X11_170412_1_2/abundances.cxb,tophat_14120X11_170414_1_2/abundances.cxb,tophat_14120X12_170412_1_2/abundances.cxb,tophat_14120X12_170414_1_2/abundances.cxb -o cuffdiff_lines/";
$output = `$cmd`;
print $output."\n";

#Cuffnorm
$cmd = "/usr/local/src/cufflinks-2.2.1.Linux_x86_64/cuffnorm -p 5 merged_asm/merged.gtf tophat_14120X1_170412_1_2/abundances.cxb,tophat_14120X1_170414_1_2/abundances.cxb,tophat_14120X2_170412_1_2/abundances.cxb,tophat_14120X2_170414_1_2/abundances.cxb,tophat_14120X3_170412_1_2/abundances.cxb,tophat_14120X3_170414_1_2/abundances.cxb,tophat_14120X4_170412_1_2/abundances.cxb,tophat_14120X4_170414_1_2/abundances.cxb,tophat_14120X5_170412_1_2/abundances.cxb,tophat_14120X5_170414_1_2/abundances.cxb,tophat_14120X6_170412_1_2/abundances.cxb,tophat_14120X6_170414_1_2/abundances.cxb tophat_14120X7_170412_1_2/abundances.cxb,tophat_14120X7_170414_1_2/abundances.cxb,tophat_14120X8_170412_1_2/abundances.cxb,tophat_14120X8_170414_1_2/abundances.cxb,tophat_14120X11_170412_1_2/abundances.cxb,tophat_14120X11_170414_1_2/abundances.cxb,tophat_14120X12_170412_1_2/abundances.cxb,tophat_14120X12_170414_1_2/abundances.cxb -o cuffnorm/";
$output = `$cmd`;
print $output."\n";


