#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

my $help = 0;
my $verbose = 1;
my $emf_directory = ".";
my $species = "homo_sapiens";


my $desc = qq{
SYNOPSIS

index_emf_files.pl --emf_directory <DIRNAME>
    
DESCRIPTION
    
This script reads all the EMF files in the directory and indexes the entries for a given
species. The output is then used by the AncestralAlleles-live plugin to locate the
alignments and run the inference on these.

OPTIONS

--emf_directory
    
    Path of the directory containing the EMF files. These can be compressed with gzip

    Default: $emf_directory

--species SPECIES
    
    The name of the species to index.

    Default: $species
    
--verbose
    
    Prints progress on the standard output.

};

GetOptions(
    "help"              => \$help,
    "verbose!"          => \$verbose,
    "emf_directory=s"   => \$emf_directory,
    "species=s"         => \$species,
);

if ($help) {
    print $desc;
    exit(0);
}

if (!-d $emf_directory) {
    die "EMF directory $emf_directory is not a directory\n";
}

my $emf_index_file = "$emf_directory/$species.index";
if (-e $emf_index_file) {
    print "Index file $emf_index_file alredy exists. Overwrite? [y/n]: ";
    my $resp = <STDIN>;
    if ($resp =~ /^\s+y/i) {
        print "Abort.";
        exit(0);
    }
}

opendir(EMFDIR, $emf_directory) or die "Cannot open EMF directory $emf_directory\n";
my @emf_files = grep {$_ =~ /\.emf(.gz|bz2)?$/i and $_ !~ /README/} readdir(EMFDIR);
closedir(EMFDIR);

my @entries =();

foreach my $this_emf_file (@emf_files) {
    print "[", scalar(localtime()), "] Parsing $this_emf_file...\n" if ($verbose);
    if ($this_emf_file =~ /\.emf\.gz$/) {
        open(EMF, "gunzip -c $emf_directory/$this_emf_file |") or die "Cannot open compressed file $this_emf_file\n";
    } elsif ($this_emf_file =~ /\.emf.bz2$/) {
        open(EMF, "bunzip2 -c $emf_directory/$this_emf_file |") or die "Cannot open compressed file $this_emf_file\n";
    } elsif ($this_emf_file =~ /\.emf$/) {
        open(EMF, "$emf_directory/$this_emf_file") or die "Cannot open file $this_emf_file\n";
    } else {
        die "I don't know how to handle a file with this extension: $this_emf_file\n";
    }

    my $emf_block_start_position = undef;
    my $current_file_position = 0; # Refers to the begining of the line
    while (<EMF>) {
        if (/^SEQ (.+)/) {
            if (!defined($emf_block_start_position)) {
                $emf_block_start_position = $current_file_position;
            } else {
                next;
            }
            my $info = $1;
            my ($species, $chromosome, $start, $end, $strand) =  $info =~ /(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/;
            if ($species eq "homo_sapiens") {
                push(@entries, [$chromosome, $start, $end, $strand, $this_emf_file, $emf_block_start_position]);
            }
        } elsif (/^DATA/) {
            $emf_block_start_position = undef;
        }
        $current_file_position = tell(EMF);
    }
    close(EMF);
}


print "[", scalar(localtime()), "] Sorting entries...\n" if ($verbose);
my @sorted_entries = sort sort_entries @entries;

print "[", scalar(localtime()), "] Writting index file...\n" if ($verbose);
open(INDEX, ">".$emf_index_file) or die "Cannot write into the index file $emf_index_file\n";
foreach my $this_entry (@sorted_entries) {
    print INDEX join("\t", @$this_entry), "\n";
}
close(INDEX);

print "[", scalar(localtime()), "] Done.\n" if ($verbose);

sub sort_entries {
    my ($chr1, $start1, $end1, $strand1, $emf1, $pos1) = @$a;
    my ($chr2, $start2, $end2, $strand2, $emf2, $pos2) = @$b;
    
    if ($chr1 eq $chr2) {
        return ($start1 <=> $start2 || $emf1 cmp $emf2 || $pos1 <=> $pos2);
    } elsif ($chr1 =~ /^\d/ and $chr2 =~ /^\d/) {
        return $chr1 <=> $chr2;
    } elsif ($chr1 =~ /^\d/) {
        return -1;
    } elsif ($chr2 =~ /^\d/) {
        return 1;
    } else {
        return $chr1 cmp $chr2;
    }
}
