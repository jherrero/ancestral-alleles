package RunAncestralAllelesOnEMF;

use strict;
use File::Spec;
use BaseAncestralAlleles;

use base ('Bio::EnsEMBL::Hive::Process');


###############################################
##
## DATA STRUCTURES
##
###############################################
##
## SORTED_ALIGNMENT
##
## This is used for the original and sub-EMF alignments. It is a hash with 3 keys:
## - tree:      a string representation of the sequence tree, with no names for the ancestral sequences
## - positions: defines the order of the sequences in the tree (required for running Ortheus)
## - sequences: a sorted array (see positions) of hashes whose keys are:
##      - species:  the species name
##      - chr:      the chromosome (or similar) name
##      - start:    the start position in that chromosome (e! coordinates)
##      - end:      the end position in that chromosomes (e! coordinates)
##      - strand:   the strand of the sequence (1 or -1)
##      - name:     a label representation of the above, as found in the tree lines of the EMF files (i.e. Hsap_1_12345_22345[+])
##      - aligned_sequence:     self-explanatory. Might not exist if the sequence is to be aligned.
##      - original_sequence:    self-explanatory. Might not exist if the sequence is already aligned.
##
###############################################
##
## ORTHEUS_ALIGNMENT
##
## This is used for capturing the Ortheus alignments. It is a hash, where each entry is one of the sequences from Ortheus.
## The key is the name of the sequence "a la Ortheus" (concatenation of the sequence numbers with an "_"), but with sorted
## leave names. Additionally, the parser also calls "ref", "sis", "anc" and "old" as aliases for the corresponding sequences.
## The values are a hash whose keys are:
## - leaves:            a hash with keys being the leaf names and values being simply 1. Useful to match leaf content.
## - num_leaves:        the number of leaves under that tree. 1 for leaves, more for ancestral nodes.
## - aligned_sequence:  self-explanatory
## - name:              the same label as for the SORTED_ALIGNMENT for leaves. Concatenation of labels these otherwise.
##
## Note: after the calling of alleles, the ortheus_alignment contain a new key-value pair with the lengths of the flank5, allele and flank3
##
###############################################


#Set of bases to insert to the left of the current base
my $inserts;
%$inserts = (
    'A' => ['A','C','G','T'],
    'C' => ['A','C','G','T'],
    'G' => ['A','C','G','T'],
    'T' => ['A','C','G','T'],
    'N' => ['A','C','G','T']);

my %event_type = %BaseAncestralAlleles::event_type;


sub param_defaults {
    
    return {
        'flank'                 => 10,
        'max_alignment_length'  => 100,
        'ortheus_exe'           => "ortheus_core",
        'nuc_frequency'         => "0.3 0.2 0.2 0.3", # Used for -n option in Ortheus (seems like not in use any longer)
        'verbose'               => 0,
        'muscle_exe'            => 0,
    }
}

sub fetch_input {
    my $self = shift @_;
    
    my $emf_file = $self->param('emf');
    my $tmp_dir = $self->worker_temp_directory;
    
    if (!-e $emf_file) {
        die "EMF file <$emf_file> does not exist\n";
    }
    
    system("cp", $emf_file, $tmp_dir);

    my ($volume, $directories, $filename) = File::Spec->splitpath($emf_file);
    # print join("\t", $volume, $directories, $filename), "\n";
    # Use the location of the copy from now on.
    $emf_file = $tmp_dir.$filename;
    
    if ($emf_file =~ /.gz$/) {
        system("gunzip", $emf_file);
        $emf_file =~ s/.gz$//;
        die "Cannot find uncompressing $emf_file" if (!-e $emf_file);
    } elsif ($emf_file =~ /.bz$/) {
        system("bunzip2", $emf_file);
        $emf_file =~ s/.bz$//;
        die "Cannot find uncompressing $emf_file" if (!-e $emf_file);
    }
    

    if ($self->param_is_defined('vcf')) {
        my $vcf_file = $self->param('vcf');
        
        my $positions;
        open(VCF, $vcf_file) or die "Cannot open vcf file $vcf_file";
        while (<VCF>) {
            my ($chr, $pos, $name, $ref_allele, $alt_allele) = split("\t", $_);
            $pos++; # VCF file includes the bp before the indel.
            $positions->{$chr}->{$pos} = 1;
        }
        close(VCF);
        $self->param("positions", $positions);
    }
    
    if ($self->param_is_defined('out')) {
        open(OUTPUT, '>', $self->param('out')) or die;
    } else {
        *OUTPUT = *STDOUT;
    }

    if ($self->param_is_defined('position')) {
        my ($chr, $pos) = split(":", $self->param('position'));
        $self->param("single_position", [$chr, $pos]);

        my $positions;
        $positions->{$chr}->{$pos} = 1;
        $self->param("positions", $positions);
    }
    
    if ($self->param_is_defined('anc')) {
        my $anc_file = $self->param('anc');
        
        my $anc_fasta;
        open(ANC, $anc_file) or die "Cannot open ancestor file $anc_file";
        while (<ANC>) {
            next if (/^>/);
            $_ =~ s/\s//g;
            $anc_fasta .= $_;
        }
        close(ANC);
        $self->param("anc_fasta", $anc_fasta);
    }
    

    $self->param("emf", $emf_file);
}

sub run {
    my $self = shift @_;
    
    my $emf_file = $self->param('emf');
    if (!-e $emf_file) {
        die "EMF file <$emf_file> does not exist\n";
    }

    open(EMF, $emf_file) or die "Cannot open EMF file $emf_file\n";
#    seek(EMF, 347826786, 0);
    my $sorted_alignment = undef;
    my $pattern = "";
    while (<EMF>) {
        if (/^SEQ (.+)/) {
            my $info = $1;
            my ($species, $chromosome, $start, $end, $strand) =  $info =~ /(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/;
            if ($species eq "ancestral_sequences") {
                # Match but don't store the sequence for the ancestral sequences
                $pattern .= " ?\\S";
            } else {
                my $name = get_name_from_data($species, $chromosome, $start, $end, $strand);

                push(@{$sorted_alignment->{sequences}}, {species=>$species, chr=>$chromosome, start=>$start, end=>$end, strand=>$strand, name=>$name});
                # Match and store the sequence for the extant sequences
                $pattern .= " ?(\\S)";
            }
        } elsif ($_ =~ /^SCORE/) {
            # Match but don't store the score
            $pattern .= " \-?[\\d\.]+";
        } elsif (/^TREE (.+)/) {
            my $tree = $1;
            $tree =~ s/\:0([^\.])/:0.0001$1/g;
            $tree =~ s/Aseq_Ancestor_[\d_]+\[.\]//g;
            $sorted_alignment->{tree} = $tree;
        } elsif (/^DATA/) {
            while(<EMF>) {
                my @this_line = uc($_) =~ /$pattern/;
                for (my $i=0; $i<@this_line; $i++) {
                    $sorted_alignment->{sequences}->[$i]->{aligned_sequence} .= $this_line[$i];
                }
                
                last if (/\/\//);
            }

            # Fix the tree.
            fix_alignment_tree($sorted_alignment);
            
            ### End of the EMF block. Send the data to the inference code.
            $self->infer_ancestral_alleles_for_this_sorted_alignment($sorted_alignment);
            
            $sorted_alignment = undef;
            $pattern="";

            ## TODO
            ## TODO
            ## TODO
            ## TODO
            ## TODO
            ## TODO
#            last;
        }
    }
}

sub write_output {
    my $self = shift @_;

    if ($self->param_is_defined('out')) {
        close(OUTPUT);
    }

}


sub fix_alignment_tree {
    my ($sorted_alignment) = @_;
    
    my $tree = $sorted_alignment->{tree};
    
    for (my $i = 0; $i < @{$sorted_alignment->{sequences}}; $i++) {
        my $name = $sorted_alignment->{sequences}->[$i]->{name};
        $tree =~ s/\Q$name\E/$i/g;
    }

    my @positions = ($tree =~ /(\d+)\:/g);
    
    $sorted_alignment->{tree} = $tree;
    
    $sorted_alignment->{positions} = \@positions;
}

sub infer_ancestral_alleles_for_this_sorted_alignment {
    my ($self, $sorted_alignment) = @_;
    
    my $work_dir = $self->worker_temp_directory;
    
    my $flank_length = $self->param('flank');
    
    my $max_alignment_length = $self->param('max_alignment_length');
    
    my $ortheus_exe = $self->param("ortheus_exe");
    
    my $muscle_exe = $self->param("muscle_exe");
    
    my $verbose = $self->param("verbose");
    
    my $positions = $self->param("positions");
    
    my $single_position;
    if ($self->param_is_defined('single_position')) {
        $single_position = $self->param("single_position");
    }
    
    my $anc_fasta;
    if ($self->param_is_defined('anc_fasta')) {
        $anc_fasta = $self->param("anc_fasta");
    }
    
    ## TODO
    ## TODO: Support for any sequence to be the reference
    ## TODO
    ## TODO
    ## TODO
    ## TODO
    ## TODO
    ## TODO

    if (!defined($sorted_alignment->{ref_seq})) {
        $sorted_alignment->{ref_seq} = 0;
    }
    my $ref_seq = $sorted_alignment->{ref_seq};

    ## Reverse-complement if required
    straighten_sorted_alignment($sorted_alignment);


    my $original_seq = $sorted_alignment->{sequences}->[$ref_seq]->{aligned_sequence};
    $original_seq =~ s/\-//g;
    my $seq_region = $sorted_alignment->{sequences}->[$ref_seq]->{chr};
    my $seq_region_start = $sorted_alignment->{sequences}->[$ref_seq]->{start};
    my $seq_region_end = $sorted_alignment->{sequences}->[$ref_seq]->{end};
    my $max_position = length($original_seq) - $flank_length * 2;

    if ($single_position) {
        print join("\t", $seq_region, $seq_region_start, $seq_region_end, $single_position->[0], $single_position->[1]), "\n";
        if ($single_position->[0] ne $seq_region or $single_position->[1] > $seq_region_end or $single_position->[1] < $seq_region_start) {
            return;
        }
    }
    
    my $pre_start = length(($sorted_alignment->{sequences}->[$ref_seq]->{aligned_sequence} =~ /^(\-*)/)[0]);
    for (my $i = 0; $i <= $max_position; $i++) {
        my $indel_coordinate = $seq_region_start + $flank_length + $i;

        my ($flank5, $flank3) = substr($sorted_alignment->{sequences}->[$ref_seq]->{aligned_sequence}, $pre_start) =~ /((?:[a-zA-Z]\-*){$flank_length})((?:\-*[a-zA-Z]){$flank_length})/;
        my $alignment_offset = $pre_start;
        my $alignment_length = length($flank5)+length($flank3);
        my $next_nucl = substr($flank3, 0, 1);

        $pre_start += length(($flank5 =~ /([a-zA-Z]\-*)/)[0]);

        ## TODO
        ## TODO: Add AA for point mutations
        ## TODO
        ## TODO
        ## TODO
        ## TODO
        ## TODO
        ## TODO
#        if ($indel_coordinate > 52520000) {
#            last;
#        }
        if ($positions) {
#            if (!defined($positions->{$seq_region}->{$indel_coordinate}) and
#                !defined($positions->{$seq_region}->{$indel_coordinate-1}) and
#                !defined($positions->{$seq_region}->{$indel_coordinate+1})
#                ) {
            if (!defined($positions->{$seq_region}->{$indel_coordinate})) {
                next;
            }
        }

        print OUTPUT join("\t", $seq_region, $indel_coordinate, $next_nucl, (substr($anc_fasta, $indel_coordinate-1, 1) or "?"), "s;");
        
#        if ($positions) {
#            if (!defined($positions->{$seq_region}->{$indel_coordinate})
#                and !defined($positions->{$seq_region}->{$indel_coordinate-1})
#                and !defined($positions->{$seq_region}->{$indel_coordinate+1})
#                ) {
#                print "\n"; # Finish the line
#                next;
#            }
#        }
        
#        print join("\t", $seq_region, $indel_coordinate, substr($anc_fasta, $indel_coordinate-1, 1), "s;");
        
        ## DEBUG: Prints information about the sub-alignment being extracted:
        ## print join(" ", $indel_coordinate, $max_position, $flank5, $flank3), "\n";

        ## DEBUG: Prints the aligned sequences of the extracted sub-alignment
        ## foreach my $this_seq (@{$sorted_alignment->{sequences}}) {
        ##     print substr($this_seq->{aligned_sequence}, $alignment_offset, $alignment_length), "\n";
        ## }

        ## DEBUG: Prints the tree and the order of the sequences
        ## print $sorted_alignment->{tree}, "\n";
        ## print join(" -- ", @{$sorted_alignment->{positions}}), "\n";

        my $reference_sequence_exception = get_reference_sequence_exception($flank5, $flank3, $alignment_length, $max_alignment_length, $work_dir, $verbose);
        if ($reference_sequence_exception) {
            print OUTPUT "e\t$reference_sequence_exception;\n";
            next;
        }
        
        
        ## Build sub-alignment for the flanking regions
        my $sub_sorted_alignment = get_sub_sorted_alignment($sorted_alignment, $alignment_offset, $alignment_length);

#        use Data::Dumper;
#        print Dumper($sub_sorted_alignment);
        splice_uninformative_sequences_from_alignment($sub_sorted_alignment);
#        print "\n\n\n";
#        print Dumper($sub_sorted_alignment);
        
        


        #####################################################################
        ##
        ## EXCEPTION: All the other sequences are all N's
        ##
        if (@{$sub_sorted_alignment->{sequences}} < 2) {
            ## No sequences left apart from the reference one.
            ## Don't mark it as an excpetion but as missing coverage only
            print OUTPUT ";e\tALL_OTHER_Ns\n";
            next;
        }
        ##
        #####################################################################
        
        
        ## DEBUG: Prints the sequences to be aligned
        ## foreach my $this_subseq (@{$sub_sorted_alignment->{sequences}}) {
        ##     print ">", $this_subseq->{name}, "\n", $this_subseq->{original_sequence}, "\n";
        ## }
        
        my $reference_ortheus_alignment = run_ortheus($sub_sorted_alignment, $ortheus_exe, $muscle_exe, $work_dir, $verbose);

        ## TODO
        ## TODO: Print the output to a file
        ## TODO
        ## TODO
        ## TODO
        ## TODO
        ## TODO
        ## TODO

        
        #####################################################################
        ##
        ## Infer the Ancestral alleles for the possible insertions (skip if same nucl as next pos, see %inserts)
        ##
        foreach my $this_indel (@{$inserts->{$next_nucl}}) {
            my $aln_seq = $flank5.$this_indel.$flank3;
            $aln_seq =~ s/\-//g;
            $sub_sorted_alignment->{sequences}->[$ref_seq]->{original_sequence} = $aln_seq;
            my $insertion_ortheus_alignment = run_ortheus($sub_sorted_alignment, $ortheus_exe, $muscle_exe, $work_dir, $verbose);
            my ($reference_allele_call, $alternate_allele_call, $ancestral_allele_call, $indel_type_call) = $self->call_ancestral_allele_for_insertion($reference_ortheus_alignment, $insertion_ortheus_alignment, $this_indel, $flank_length, $verbose);
            print OUTPUT join("\t", $this_indel, "i", ($event_type{$indel_type_call} or $indel_type_call), $reference_allele_call, $alternate_allele_call, $ancestral_allele_call.";");
        }
        ##
        #####################################################################
        
        
        #####################################################################
        ##
        ## Infer the Ancestral alleles for the deletion
        ##
        my $aln_seq = $flank5.substr($flank3, 1);
        $aln_seq =~ s/\-//g;
        $sub_sorted_alignment->{sequences}->[$ref_seq]->{original_sequence} = $aln_seq;
        my $deletion_ortheus_alignment = run_ortheus($sub_sorted_alignment, $ortheus_exe, $muscle_exe, $work_dir, $verbose);
        my ($reference_allele_call, $alternate_allele_call, $ancestral_allele_call, $indel_type_call) = $self->call_ancestral_allele_for_deletion($reference_ortheus_alignment, $deletion_ortheus_alignment, $next_nucl, $flank_length, $verbose);
        print OUTPUT join("\t", $next_nucl, "d", ($event_type{$indel_type_call} or $indel_type_call), $reference_allele_call, $alternate_allele_call, $ancestral_allele_call.";");
        ##
        #####################################################################
        
        #####################################################################
        ##
        ## Finish line
        ##
        print OUTPUT "\n";
        ##
        #####################################################################
#        last;
#        last if ($i > 1000);
    }
}


1;