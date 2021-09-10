#!/usr/bin/perl -w

use warnings;
use strict;
use POSIX;
use Getopt::Long;
use File::Path;
use File::Basename;
use Pod::Usage;

=head1 NAME

Generate bfGWAS Makefile

=head1 SYNOPSIS

./gen_mkf.pl [options] 

Options:
  -h      help messages
  -d      degub
  --mf    output make file

=head1 OPTIONS

=over 8

=item B<-help>

  --wkdir         work directory : location for all output files
  --tool_dir         directory for the C++ executable file for Estep (running MCMC)
  --LDdir       directory of LD correlation files
  --Score_dir    directory of GWAS single variant score statistic files
  --filehead         text file with a list of fileheads for genome blocks
  --anno_dir       annotation file directory
  --anno_code       annotation classification code file
  --hfile      initial hyper parameter value file
  --pv       phenotype variance
  --Nsample       sample size
  --maf      maf threshold: default 0.5% 
  --Nburnin         number of burn ins: default 10,000
  --Nmcmc         number of MCMC iterations: default 10,000
  --NmcmcLast       number of MCMC iterations for the last EM iteration: default 50,000
  --pp       specify prior for the causal probability: default 1e-6
  --abgamma  specify inverse gamma prior for the effect size variance: default 0.1
  --win      window size of the neighborhood: default 100
  -c         compress genotype data (1) or not (0): default 0 (do not compress)
  --initype  specify initial model (1:start with top signal), (2: start with genome-wide significant signals), or default (3: stepwise selected signals)
  --smin     minimum number of variates per block in the model: default 0
  --smax     maximum number of variates per block in the model: default 5
  --em       number of EM iterations: default 3
  --mf       output make file

=back

=head1 DESCRIPTION

B<gen_mkf.pl> will generate a makefile for conducting Bayesian Functional GWASs by bfGWAS. 

=cut

# define default option variables
my $help;
my $verbose;
my $debug;
my $man;
my $launchMethod = "local";
my $makeFile = "BFGWAS.mk";

my $tool_dir="";
my $wkdir=getcwd();
my $filehead = "";
my $Score_dir = "";
my $LDdir="";
my $Nsample;
my $anno_dir="";
my $annoCode="";
my $hfile="";

my $EM=3;
my $maf="0.005";
my $smin="0";
my $smax="5";
my $win="100";
my $burnin="10000";
my $Nmcmc="10000";
my $NmcmcLast="50000";
my $compress=0;
my $initype="3";
my $pv="1";
my $pp="1e-6";
my $abgamma="0.1";

#initialize options
Getopt::Long::Configure ('bundling');

if(!GetOptions ('h'=>\$help, 'v'=>\$verbose, 'd'=>\$debug, 'm'=>\$man,
                'wkdir:s'=>\$wkdir, 'tool_dir:s' =>\$tool_dir, 
                'LDdir:s'=>\$LDdir, 'Score_dir:s'=>\$Score_dir, 
                'filehead:s'=>\$filehead, 'anno_dir:s'=>\$anno_dir, 
                'anno_code:s'=>\$annoCode, 'hfile:s'=>\$hfile, 
                'pv:s'=>\$pv,'Nsample:i'=>\$Nsample, 'maf:s'=>\$maf,
                'Nburnin:s'=>\$burnin, 'Nmcmc:s'=>\$Nmcmc, 'NmcmcLast:s'=>\$NmcmcLast,
                'pp:s'=>\$pp, 'abgamma:s'=>\$abgamma, 'win:s'=>\$win,
                'c:i'=>\$compress, 'initype:s'=>\$initype,
                'smin:s'=>\$smin, 'smax:s'=>\$smax, 
                'em:i'=>\$EM, 'mf:s'=>\$makeFile)
  || !defined($wkdir) || scalar(@ARGV)!=0)
{
    if ($help)
    {
        pod2usage(1);
        exit(0);
    }
    elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }
    else
    {
        pod2usage(1);
        exit(0);
    }
}

if ($help)
    {
        pod2usage(1);
        exit(0);
    }
elsif($man) {
        pod2usage(-exitval=>0, -verbose => 2);
    }


my $toolE="${tool_dir}/bin/Estep_mcmc";
my $rs="${tool_dir}/bin/Mstep.r";

##############
#print options
##############
printf("Options\n");
printf("\n");
printf("launch method : %s\n", $launchMethod);
printf("work directory : %s\n", $wkdir);
print "Estep: ", $toolE, "\n", "Rscript: ", $rs, "\n"; 
print "Score_dir: ", $Score_dir, "\n", "LDdir: ", $LDdir, "\n", 
        "anno_dir: ", $anno_dir, "\nannoCode: ", $annoCode, "\n", 
        "hfile: ", $hfile, "\nfilehead text file: ", $filehead, "\n",
        "maf ", $maf, "; smin ", $smin, "\n", 
        "smax ", $smax, "; win ", $win, "\n", 
        "burnin ", $burnin, "; Nmcmc ", $Nmcmc, "; NmcmcLast ", $NmcmcLast, "\n",
        "compress ", $compress, "; initype ", $initype, "\n",
        "pv ", $pv, "; pp ", $pp, "; abgamma ", $abgamma, "\n";
printf("\n");

my $comp = "";
if ($compress != 0){
  $comp = "-comp"
}

#arrays for storing targets, dependencies and commands
my @tgts = ();
my @deps = ();
my @cmds = ();

#temporary variables
my $tgt;
my $dep;
my @cmd;

mkpath($wkdir);

########################################
# Initial Set up before EM iterations
########################################

### prepare files before EM_MCMC
my $hypcurrent="$wkdir/hypval.current";
$tgt = "$wkdir/pre_em.OK";
$dep = "";
@cmd = "rm -f -r $wkdir/output $wkdir/Eoutput $wkdir/OUT";
push(@cmd, "mkdir -p $wkdir/output $wkdir/Eoutput $wkdir/OUT");
push(@cmd, "cp -f $hfile $hypcurrent");
push(@cmd, "> $wkdir/Eoutput/EM\_result.txt");
push(@cmd, "> $wkdir/Rout.txt");
makeJob($tgt, $dep, @cmd);  


###### EM step 0 without dependencies ###########
my $i;
my $premcmcOK="";
my @filehead;
my $line;

open(my $FILELIST, $filehead)
    or die "Can not open $filehead \!";
 while ($line = <$FILELIST>) {
    chop $line;
    push(@filehead, $line);
}
close $FILELIST;

if(@filehead == 0) {
    print STDERR "file list is empty\! \n";
    exit(1);
} 
else{ print "Total \# of genome blocks: ", scalar(@filehead), "\n \n"; }


for(my $j=0; $j< @filehead; ++$j)
    {
        $i=0;
        $line = $filehead[$j];
        $premcmcOK .= "$wkdir/OUT/$line.$i.OK ";
        $tgt = "$wkdir/OUT/$line.$i.OK";
        $dep = "$wkdir/pre_em.OK";

        @cmd = "$toolE -inputSS -score ${Score_dir}/${line}.score.txt.gz -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a ${anno_dir}/Anno_${line}.txt.gz -fcode ${annoCode} -hfile ${hypcurrent} -n ${Nsample} -pv ${pv} -maf ${maf} -bvsrm -smin $smin -smax $smax -win $win -o ${line} -w ${burnin} -s ${Nmcmc} -initype $initype -seed 2021 > ${wkdir}/OUT/${line}.output.txt";
        makeJob($tgt, $dep, @cmd);
    }

my $paramfile="$wkdir/Eoutput/paramtemp$i.txt";
my $hypfile="$wkdir/Eoutput/hyptemp$i.txt";
my $logfile="$wkdir/Eoutput/log$i.txt";

$tgt = "$wkdir/Eoutput/cp_param$i.OK";
$dep = "$premcmcOK";
@cmd = "cat \`ls -d -1 $wkdir/output/** | grep paramtemp | sort\` > $paramfile";
push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep hyptemp | sort\` > $hypfile");
push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep log | sort\` > $logfile");
makeJob($tgt, $dep, @cmd);

$tgt = "$wkdir/R$i.OK";
$dep = "$wkdir/Eoutput/cp_param$i.OK $wkdir/pre_em.OK";
@cmd = "Rscript --vanilla ${rs} $hypfile $i $pp $abgamma $wkdir/Eoutput/EM_result.txt $hypcurrent >> $wkdir/Rout.txt";
makeJob($tgt, $dep, @cmd);


####### With dependencies of previous output 
my $ipre="";

for $i (1..$EM){

    $ipre=$i-1; $premcmcOK="";

    for(my $j=0; $j< @filehead; ++$j){
        $line=$filehead[$j];
        $premcmcOK .= "$wkdir/OUT/$line.$i.OK ";
        $tgt = "$wkdir/OUT/$line.$i.OK";
        $dep = "$wkdir/R$ipre.OK";
        if($i < $EM){
          @cmd = "$toolE -inputSS -score ${Score_dir}/${line}.score.txt.gz -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a ${anno_dir}/Anno_${line}.txt.gz -fcode ${annoCode} -hfile ${hypcurrent} -n ${Nsample} -pv ${pv} -maf ${maf} -bvsrm -smin $smin -smax $smax -win $win -o ${line} -w ${burnin} -s ${Nmcmc} -initype $initype -seed 2021 > ${wkdir}/OUT/${line}.output.txt" ;
        } elsif ($i == $EM){
            @cmd = "$toolE -inputSS -score ${Score_dir}/${line}.score.txt.gz -LDcorr ${LDdir}/${line}.LDcorr.txt.gz -a ${anno_dir}/Anno_${line}.txt.gz -fcode ${annoCode} -hfile ${hypcurrent} -n ${Nsample} -pv ${pv} -maf ${maf} -bvsrm -smin $smin -smax $smax -win $win -o ${line} -w ${burnin} -s ${NmcmcLast} -initype $initype -seed 2021 > ${wkdir}/OUT/${line}.output.txt";
      }
        makeJob($tgt, $dep, @cmd);
    }

    $paramfile="$wkdir/Eoutput/paramtemp$i.txt";
    $hypfile="$wkdir/Eoutput/hyptemp$i.txt";
    $logfile="$wkdir/Eoutput/log$i.txt";

  $tgt = "$wkdir/Eoutput/cp_param$i.OK";
  $dep = "$premcmcOK";
  @cmd = "cat \`ls -d -1 $wkdir/output/** | grep paramtemp | sort \` > $paramfile";
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep hyptemp | sort\` > $hypfile");
  push(@cmd, "cat \`ls -d -1 $wkdir/output/** | grep log | sort\` > $logfile");
  makeJob($tgt, $dep, @cmd);

  $tgt = "$wkdir/R$i.OK";
  $dep = "$wkdir/Eoutput/cp_param$i.OK";
  @cmd = "Rscript --vanilla $rs $hypfile $i $pp $abgamma $wkdir/Eoutput/EM_result.txt $hypcurrent >> $wkdir/Rout.txt";
  makeJob($tgt, $dep, @cmd);

}


#*******************
#Write out make file
#*******************
open(MAK,">$makeFile") || die "Cannot open $makeFile\n";
#this tells GNU Make to delete files that were generated during a command that did not end successfully.
print MAK ".DELETE_ON_ERROR:\n\n";
#this ensures that all the targets are executed; exclude clean
print MAK "all: @tgts\n\n";

######## Create clean jobs command #######
push(@tgts, "clean");
push(@deps, "");
push(@cmds, "\t-rm -rf $wkdir/*.OK $wkdir/Eoutput/*.OK $wkdir/OUT/*.OK $wkdir/slurm_err/*.err");

for(my $i=0; $i < @tgts; ++$i)
{
    print MAK "$tgts[$i]: $deps[$i]\n";
    print MAK "$cmds[$i]\n";
}
close MAK;

##########
#functions
##########

#run a job locally
sub makeJob
{    
    my ($tgt, $dep, @cmd) = @_;

    push(@tgts, $tgt);
    push(@deps, $dep);
    my $cmd = "";
    for my $c (@cmd)
    {
        $cmd .= "\t" . $c . "\n";
    }
    $cmd .= "\ttouch $tgt\n";
    push(@cmds, $cmd);
}


# Check empty directories
sub dir_is_empty
{
    my ($path) = @_;
    opendir DIR, $path;
    while(my $entry = readdir DIR) {
        next if($entry =~ /^\.\.?$/);
        closedir DIR;
        return 0;
    }
    closedir DIR;
    return 1;
}
