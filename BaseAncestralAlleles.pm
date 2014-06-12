package BaseAncestralAlleles;

use strict;
use warnings;

=pod

=head1 MODULE

BaseAncestralAlleles

=head1 DESCRIPTION

This module contain the internal methods used for polarising indels in a genome. It is meant to be used as a BASE class for
the RunAncestralAlleles* modules. These will handle the input and output layer.

=head1 Exported methods

This module exports by default the following methods:

=over

=item straighten_sorted_alignment

=item get_name_from_data

=item get_sub_sorted_alignment

=item splice_uninformative_sequences_from_alignment

=item get_reference_sequence_exception

=item run_ortheus

=item call_ancestral_allele_for_deletion

=item call_ancestral_allele_for_insertion

=back

=head1 INTERNAL DATA STRUCTURES

=head2 SORTED_ALIGNMENT

 This is used for the original and sub-EMF alignments. It is a hash with 3 keys:
 - tree:      a string representation of the sequence tree, with no names for the ancestral sequences
 - positions: defines the order of the sequences in the tree (required for running Ortheus)
 - sequences: a sorted array (see positions) of hashes whose keys are:
      - species:  the species name
      - chr:      the chromosome (or similar) name
      - start:    the start position in that chromosome (e! coordinates)
      - end:      the end position in that chromosomes (e! coordinates)
      - strand:   the strand of the sequence (1 or -1)
      - name:     a label representation of the above, as found in the tree lines of the EMF files (i.e. Hsap_1_12345_22345[+])
      - aligned_sequence:     self-explanatory. Might not exist if the sequence is to be aligned.
      - original_sequence:    self-explanatory. Might not exist if the sequence is already aligned.

=head2 ORTHEUS_ALIGNMENT

 This is used for capturing the Ortheus alignments. It is a hash, where each entry is one of the sequences from Ortheus.
 The key is the name of the sequence "a la Ortheus" (concatenation of the sequence numbers with an "_"), but with sorted
 leave names. Additionally, the parser also calls "ref", "sis", "anc" and "old" as aliases for the corresponding sequences.
 The values are a hash whose keys are:
 - leaves:            a hash with keys being the leaf names and values being simply 1. Useful to match leaf content.
 - num_leaves:        the number of leaves under that tree. 1 for leaves, more for ancestral nodes.
 - aligned_sequence:  self-explanatory
 - name:              the same label as for the SORTED_ALIGNMENT for leaves. Concatenation of labels these otherwise.

 Note: after the calling of alleles, the ortheus_alignment contain a new key-value pair with the lengths of the flank5, allele and flank3

=head2 TREES

 The trees are represented as a nested set of tree nodes. The keys are:
 - children: an arrayref to a set of other tree nodes
 - distance_to_parent: a number representing the distance to the parent node in the tree
 - label: a string representing the name of this node (only for leaves)
 - mark: additional flag used when trimming the trees

=head1 INTERNAL METHODS

The rest of the documentation refers to the internal methods implemented in this module.

=cut


use Exporter 'import';

our @EXPORT = qw(
straighten_sorted_alignment
get_name_from_data
get_sub_sorted_alignment
splice_uninformative_sequences_from_alignment
get_reference_sequence_exception
run_ortheus
call_ancestral_allele_for_deletion
call_ancestral_allele_for_insertion
);

our %event_type = (
    'deletion' => 1,
    'complex_deletion' => 2,
    'funny_deletion' => 3,
    'insertion' => 4,
    'complex_insertion' => 5,
    'funny_insertion' => 6,
    'unsure' => 7);


sub get_name_from_data {
    my ($species, $chromosome, $start, $end, $strand) = @_;
    
    my $name = $species;
    $name =~ s/(.).+_([^_]{3})[^_]*/\u$1$2/;
    $name = join("_", $name, $chromosome, $start, $end)."[".($strand==-1?"-":"+")."]";
    
    return $name;
}



sub get_reference_sequence_exception {
    my ($flank5, $flank3, $alignment_length, $max_alignment_length, $work_dir, $verbose) = @_;
   
    ## EXCEPTION: the alignment is too long
    if ($alignment_length > $max_alignment_length) {
        return "LONG_INSERTION";
    }
    
    ## EXCEPTION: Either flank is Ns only
    if ($flank5 =~ /^[Nn\-]+$/ or $flank3 =~ /^[Nn\-]+$/) {
        return "ALL_N";
    }
    
#    return undef;

    ## EXCEPTION: Low-complexity sequence. Skip if either flank contains
    ## only 1 different nucleotide (not counting Ns)
    if ($flank5 =~ /[Aa]/ + $flank5 =~ /[Cc]/ + $flank5 =~ /[Gg]/ + $flank5 =~ /[Tt]/ < 2) {
        return "LOW_COMPLEXITY (HR)";
    }
    if ($flank3 =~ /[Aa]/ + $flank3 =~ /[Cc]/ + $flank3 =~ /[Gg]/ + $flank3 =~ /[Tt]/ < 2) {
        return "LOW_COMPLEXITY (HR)";
    }

    my $exception = check_complexity_of_sequence($flank5.$flank3, $work_dir);
    return $exception if ($exception);

    return undef;
}


sub check_complexity_of_sequence {
    my ($flank, $work_dir) = @_;
    
    $flank =~ s/\-//g;
#    $flank = substr($flank, 1, -1);
    $flank = substr($flank, 2, -2);
#    $flank = substr($flank, 3, -3);

    my $muscle_in = "$work_dir/flank.fa";
    
    for (my $offset = 2; $offset <= 4; $offset++) {
        open(MUSCLE, ">$muscle_in") or die();
        print MUSCLE ">locus1.$offset\n", substr($flank, $offset), "\n";
        print MUSCLE ">locus2.$offset\n", substr($flank, 0, -$offset), "\n";
        close(MUSCLE);
        my $muscle_cmd = "muscle -maxiters 10 -diags -in $muscle_in -quiet";
        my @muscle_lines = qx"$muscle_cmd";
        if ($muscle_lines[1] ne ("-"x$offset).substr($flank, $offset)."\n" or
            $muscle_lines[3] ne substr($flank, 0, -$offset).("-"x$offset)."\n") {
                my $num_matches = 0;
                for (my $pos = 0; $pos < length($muscle_lines[1])-1; $pos++) {
                    $num_matches++ if (substr($muscle_lines[1], $pos, 1) eq substr($muscle_lines[3], $pos, 1));
                }
                if ($num_matches > length($flank) - $offset - 3) {
#                if ($num_matches > length($flank) - $offset - 5) {
#                    print "\n", @muscle_lines,"matches:$num_matches\n\n";
                    unlink($muscle_in);
                    return "LOW_COMPLEXITY (muscle $offset $num_matches ". (length($flank) - $offset - 3). " $flank)";
                }
            }
    }
    unlink($muscle_in);
    return undef;
}

sub run_ortheus {
    my ($sorted_alignment, $ortheus_exe, $muscle_exe, $work_dir, $verbose) = @_;
    my $ortheus_alignment;
    
    my $tree = $sorted_alignment->{tree};
    my $fasta_filenames = [];
    my $muscle_in = "$work_dir/muscle.in.$$.fa";
    my $muscle_out = "$work_dir/muscle.out.$$.fa";
    my $muscle_order;
    
    if ($muscle_exe) {
        open(MUSCLE, ">$muscle_in") or die;
    }
    
    for (my $i = 0; $i < @{$sorted_alignment->{positions}}; $i++) {
        my $this_seq = $sorted_alignment->{sequences}->[$sorted_alignment->{positions}->[$i]];
        my $file_name = "$work_dir/ortheus.$$.".$sorted_alignment->{positions}->[$i].".fa";
        open(FASTA, ">$file_name") or die "Cannot open $file_name";
        print FASTA ">", $this_seq->{name}, "\n", $this_seq->{original_sequence}, "\n";
        if ($muscle_exe) {
            print MUSCLE ">", $this_seq->{name}, "\n", $this_seq->{original_sequence}, "\n";
            $muscle_order->{$this_seq->{name}} = $i;
        }
        close(FASTA);
        push(@$fasta_filenames, $file_name);
    }
    if ($muscle_exe) {
        close(MUSCLE);
        
        my $muscle_cmd = "muscle -maxiters 1 -diags -in $muscle_in -out $muscle_out";
        system($muscle_cmd);
        my @muscle_lines;
        open(MUSCLE, $muscle_out);
        while (<MUSCLE>) {
            if (/>(.+)/) {
                my $order = $muscle_order->{$1};
                $muscle_lines[$order*2] = $_;
                $muscle_lines[$order*2+1] = <MUSCLE>;
            }
        }
        close(MUSCLE);
        open(MUSCLE, ">".$muscle_out);
        print MUSCLE join("", @muscle_lines);
        close(MUSCLE);
    }
    
    
    my $output_alignment_file = "$work_dir/ortheus.$$.output.mfa";
    my $output_score_file = "$work_dir/ortheus.$$.score.txt";
    
    my $ortheus_cmd = "$ortheus_exe -b '$tree' -a ".join(" ", @$fasta_filenames)." -i 10 -j 0 -d $output_alignment_file -x $output_score_file";
    
    my @ortheus_args = ("-b", $tree, "-a", @$fasta_filenames, "-d", $output_alignment_file, "-x", $output_score_file);
    if ($muscle_exe) {
        push(@ortheus_args, "-c", $muscle_out);
        $ortheus_cmd .= " -c $muscle_out";
    }
#    system($ortheus_exe, @ortheus_args) == 0 or die "system call to ortheus failed (".join(" ", @ortheus_args).": $?";
    qx"$ortheus_cmd";
    
    open(MFA, $output_alignment_file) or die "ORTHEUS: $ortheus_cmd\nCannot open <$output_alignment_file>";
    my $sort_code;
    while(<MFA>) {
        if (/^>(\S+)/) {
            my $seq_name = $1;
            my @leaves = split("_", $seq_name);
            $sort_code = join("_", sort {$a <=> $b} @leaves);
            foreach my $leaf (@leaves) {
                $ortheus_alignment->{$sort_code}->{leaves}->{$leaf} = 1;
            }
            $ortheus_alignment->{$sort_code}->{num_leaves} = @leaves;
            if (@leaves == 1) {
                $ortheus_alignment->{$sort_code}->{name} = $sorted_alignment->{sequences}->[$sort_code]->{name};
            }
        } elsif (defined($sort_code)) {
            chomp($_);
            $ortheus_alignment->{$sort_code}->{aligned_sequence} .= $_;
        }
    }
    close(MFA);
    
    ## Find the ref, anc, old and sis sequences in the ortheus tree
    my ($ref_seq_name, $anc_seq_name, $old_seq_name, $sis_seq_name);
    # Sort the nodes by number of leaves to get the reference, then the ancestral and finally the older sequence
    foreach my $seq_name (sort {$ortheus_alignment->{$a}->{num_leaves} <=> $ortheus_alignment->{$b}->{num_leaves}} keys %$ortheus_alignment) {
        # Skip if it doesn't include the reference sequence
        next if (!$ortheus_alignment->{$seq_name}->{leaves}->{0});
        if (!defined($ref_seq_name)) {
            $ref_seq_name = $seq_name;
        } elsif (!defined($anc_seq_name)) {
            $anc_seq_name = $seq_name;
            # The code name for the sister sequence must be the same as the ancestral minus the initial "0_"
            $sis_seq_name = $seq_name;
            $sis_seq_name =~ s/^0_//;
        } elsif (!defined($old_seq_name)) {
            $old_seq_name = $seq_name;
            last; # No point in going any further
        }
    }
    
    ## It is OK not to have an "older" sequence. This can happen for shallow alignment or for cases where the reference is an out-sequence.
    if (!defined($ref_seq_name) or !defined($sis_seq_name) or !defined($anc_seq_name)) {
        print_ortheus_alignment($ortheus_alignment);
        die "Problem parsing Ortheus output (".join(" ", @ortheus_args). ")";
    }
    
    $ortheus_alignment->{"ref"} = $ortheus_alignment->{$ref_seq_name};
    $ortheus_alignment->{"sis"} = $ortheus_alignment->{$sis_seq_name};
    $ortheus_alignment->{"anc"} = $ortheus_alignment->{$anc_seq_name};
    $ortheus_alignment->{"old"} = $ortheus_alignment->{$old_seq_name} if ($old_seq_name);
    
    #    print_ortheus_alignment($ortheus_alignment);
    
    return $ortheus_alignment;
}


sub call_ancestral_allele_for_deletion {
    my ($self, $reference_ortheus_alignment, $deletion_ortheus_alignment, $this_deletion, $flank_length, $verbose) = @_;
    
    my ($ref_reference_allele, $sis_reference_allele, $anc_reference_allele, $old_reference_allele) = get_alleles_for_deletion($reference_ortheus_alignment, $flank_length, $this_deletion);
    my ($ref_deletion_allele, $sis_deletion_allele, $anc_deletion_allele, $old_deletion_allele) = get_alleles_for_deletion($deletion_ortheus_alignment, $flank_length, $this_deletion);
    
    ## TODO
    ## TODO: Add confidence information
    ## TODO
    ## TODO
    ## TODO
    ## TODO
    ## TODO
    
    my $ancestral_allele_call = $anc_reference_allele;
    my $indel_call;
    if (uc($anc_reference_allele) eq uc($anc_deletion_allele)) {
        if (uc($anc_reference_allele) eq uc($ref_reference_allele)) {
            $indel_call = "deletion";
        } elsif (uc($anc_reference_allele) eq uc($ref_deletion_allele)) {
            $indel_call = "insertion";
        } elsif (length($anc_reference_allele) == length($ref_reference_allele)) {
            $indel_call = "complex_deletion";
        } elsif (length($anc_reference_allele) == length($ref_deletion_allele)) {
            $indel_call = "complex_insertion";
        } else {
            $indel_call = "cryptic_indel";
        }
    } elsif (length($anc_reference_allele) eq length($anc_deletion_allele)) {
        if (uc($anc_reference_allele) eq uc($ref_reference_allele)) {
            $indel_call = "deletion";
        } elsif (uc($anc_deletion_allele) eq uc($ref_deletion_allele)) {
            $ancestral_allele_call = $anc_deletion_allele;
            $indel_call = "insertion";
        } elsif (length($anc_reference_allele) == length($ref_reference_allele)) {
            $indel_call = "complex_deletion";
        } elsif (length($anc_deletion_allele) == length($ref_deletion_allele)) {
            $ancestral_allele_call = $anc_deletion_allele;
            $indel_call = "complex_insertion";
        } else {
            $indel_call = "cryptic_indel";
        }
    } else {
        $ancestral_allele_call = "?";
        $indel_call = "unsure";
    }
    
    $ref_reference_allele ||= "-";
    $ref_deletion_allele ||= "-";
    $ancestral_allele_call ||= "-";
    
    if ($verbose) {
        print "\nDELETION (-$this_deletion)\n";
        print "Alignment for the reference allele:\n";
        print_ortheus_alignment($reference_ortheus_alignment);
        print "Alignment for the alternate allele (-$this_deletion):\n";
        print_ortheus_alignment($deletion_ortheus_alignment);
        print "ANCESTRAL ALLELE <$ancestral_allele_call>\n";
        print "TYPE: $indel_call\n";
    }
    
    return (($ref_reference_allele or "-"), ($ref_deletion_allele or "-"), ($ancestral_allele_call or "-"), $indel_call);
}


sub call_ancestral_allele_for_insertion {
    my ($self, $reference_ortheus_alignment, $insertion_ortheus_alignment, $this_insertion, $flank_length, $verbose) = @_;
    
    my ($ref_reference_allele, $sis_reference_allele, $anc_reference_allele, $old_reference_allele) = get_alleles_for_insertion($reference_ortheus_alignment, $flank_length, $this_insertion);
    my ($ref_insertion_allele, $sis_insertion_allele, $anc_insertion_allele, $old_insertion_allele) = get_alleles_for_insertion($insertion_ortheus_alignment, $flank_length, $this_insertion);
    
    
    ## TODO
    ## TODO: Add confidence information
    ## TODO
    ## TODO
    ## TODO
    ## TODO
    ## TODO
    
    my $ancestral_allele_call = $anc_reference_allele;
    my $indel_call;
    if (uc($anc_reference_allele) eq uc($anc_insertion_allele)) {
        if (uc($anc_reference_allele) eq uc($ref_reference_allele)) {
            $indel_call = "insertion";
        } elsif (uc($anc_reference_allele) eq uc($ref_insertion_allele)) {
            $indel_call = "deletion";
        } elsif (length($anc_reference_allele) == length($ref_reference_allele)) {
            $indel_call = "complex_insertion";
        } elsif (length($anc_reference_allele) == length($ref_insertion_allele)) {
            $indel_call = "complex_deletion";
        } else {
            $indel_call = "cryptic_indel";
        }
    } elsif (length($anc_reference_allele) eq length($anc_insertion_allele)) {
        if (uc($anc_reference_allele) eq uc($ref_reference_allele)) {
            $indel_call = "insertion";
        } elsif (uc($anc_insertion_allele) eq uc($ref_insertion_allele)) {
            $ancestral_allele_call = $anc_insertion_allele;
            $indel_call = "deletion";
        } elsif (length($anc_reference_allele) == length($ref_reference_allele)) {
            $indel_call = "complex_insertion";
        } elsif (length($anc_reference_allele) == length($ref_insertion_allele)) {
            $ancestral_allele_call = $anc_insertion_allele;
            $indel_call = "complex_deletion";
        } else {
            $indel_call = "cryptic_indel";
        }
    } else {
        $ancestral_allele_call = "?";
        $indel_call = "unsure";
    }
    if ($verbose) {
        print "\nINSERTION ($this_insertion)\n";
        print "Alignment for the reference allele:\n";
        print_ortheus_alignment($reference_ortheus_alignment);
        print "Alignment for the alternate allele ($this_insertion):\n";
        print_ortheus_alignment($insertion_ortheus_alignment);
        print "ANCESTRAL ALLELE <$ancestral_allele_call>\n";
        print "TYPE: $indel_call\n";
    }
    
    return (($ref_reference_allele or "-"), ($ref_insertion_allele or "-"), ($ancestral_allele_call or "-"), $indel_call);
}


=head2 get_alleles_for_deletion
 
 Arg[1]:        hashref $ortheus_alignment
 Arg[2]:        int $flank_length
 Arg[3]:        char $indel
 Returns:       array with string $ref_allele, string $sis_allele, string $anc_allele, string
                $old_allele, int $length_flank5, int $length_allele, int $length_flank3.
 Description:   Knowing the $flank_length and the $indel, this method parses the $ortheus_alignment
                to extract the ref, sis, anc and old alleles. It also returns the lengths of the
                flanks and the allele. The allele might be longer than 1bp is this deletion affects
                a homopolymer run.
 
                The method assumes that the $ortheus_alignment already knows the "ref", "sis", "anc"
                and "old" sequences (this is usually done by the run_ortheus() method).
 
                Note that this method works for both the reference and the alternate allele ortheus
                alignment. The expectation for most cases is to get:
                REF.ALLELE.ALIGN => ("A", "A", "A", "A", 10, 1, 9);
                ALT.ALLELE.ALIGN => ("", "A", "A", "A", 10, 1, 9);
 
=cut

sub get_alleles_for_deletion {
    my ($ortheus_alignment, $flank_length, $indel) = @_;
    
    
    my $flank3_length = $flank_length - 1;
    my ($flank5, $flank3) = $ortheus_alignment->{"ref"}->{aligned_sequence} =~ /((?:\-*[a-zA-Z]){$flank_length}).*((?:[a-zA-Z]\-*){$flank3_length})$/;

    if (!$flank5 or !$flank3) {
        print_ortheus_alignment($ortheus_alignment);
        die "Cannot extract the flanks from sequence (deletion): ".$ortheus_alignment->{"ref"}->{aligned_sequence};
    }
    
    $flank5 =~ s/\-*($indel\-*)+$//g;
    my $length_flank5 = length($flank5);
    $flank3 =~ s/^(\-*$indel)+\-*//g;
    my $length_flank3 = length($flank3);
    my $length_allele = length($ortheus_alignment->{"ref"}->{aligned_sequence}) - $length_flank5 - $length_flank3;
    $ortheus_alignment->{"lengths"} = [$length_flank5, $length_allele, $length_flank3];
    my $ref_allele = substr($ortheus_alignment->{"ref"}->{aligned_sequence}, $length_flank5, $length_allele);
    my $sis_allele = substr($ortheus_alignment->{"sis"}->{aligned_sequence}, $length_flank5, $length_allele);
    my $anc_allele = substr($ortheus_alignment->{"anc"}->{aligned_sequence}, $length_flank5, $length_allele);
    my $old_allele;
    if ($ortheus_alignment->{"old"}) {
        $old_allele = substr($ortheus_alignment->{"old"}->{aligned_sequence}, $length_flank5, $length_allele);
    }
    #    print_ortheus_alignment($ortheus_alignment);
    
    for (my $pos = 0; $pos < length($anc_allele); $pos++) {
        my $anc = substr($anc_allele, $pos, 1);
        my $sis = substr($sis_allele, $pos, 1);
        my $old = "";
        if ($old_allele) {
            $old = substr($old_allele, $pos, 1);
        }
        if ($anc eq $sis and $anc eq $old) {
            $anc = uc($anc)
        } elsif ($anc eq $sis or $anc eq $old) {
            $anc = lc($anc);
        } else {
            $anc = "?";
        }
        substr($anc_allele, $pos, 1, $anc);
    }
    
    $ref_allele =~ s/\-//g;
    $sis_allele =~ s/\-//g;
    $anc_allele =~ s/\-//g;
    $old_allele =~ s/\-//g if ($old_allele);
    
    return($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3);
}


=head2 get_alleles_for_insertion
 
 Arg[1]:        hashref $ortheus_alignment
 Arg[2]:        int $flank_length
 Arg[3]:        char $indel
 Returns:       array with string $ref_allele, string $sis_allele, string $anc_allele, string
                $old_allele, int $length_flank5, int $length_allele, int $length_flank3.
 Description:   Knowing the $flank_length and the $indel, this method parses the $ortheus_alignment
                to extract the ref, sis, anc and old alleles. It also returns the lengths of the
                flanks and the allele. The allele might be longer than 1bp is this insertion affects
                a homopolymer run.
 
                The method assumes that the $ortheus_alignment already knows the "ref", "sis", "anc"
                and "old" sequences (this is usually done by the run_ortheus() method).
 
                Note that this method works for both the reference and the alternate allele ortheus
                alignment. The expectation for most cases is to get:
                REF.ALLELE.ALIGN => ("", "", "", "", 10, 0, 10);
                ALT.ALLELE.ALIGN => ("A", "", "", "", 10, 1, 10);

=cut

sub get_alleles_for_insertion {
    my ($ortheus_alignment, $flank_length, $indel) = @_;
    
    my ($flank5, $flank3) = $ortheus_alignment->{"ref"}->{aligned_sequence} =~ /((?:\-*[a-zA-Z]){$flank_length}).*((?:[a-zA-Z]\-*){$flank_length})$/;
    if (!$flank5 or !$flank3) {
        print_ortheus_alignment($ortheus_alignment);
        die "Cannot extract the flanks from sequence (insertion): ".$ortheus_alignment->{"ref"}->{aligned_sequence};
    }
    
    $flank5 =~ s/\-*($indel\-*)+$//g;
    my $length_flank5 = length($flank5);
    $flank3 =~ s/^(\-*$indel)+\-*//g;
    my $length_flank3 = length($flank3);
    my $length_allele = length($ortheus_alignment->{"ref"}->{aligned_sequence}) - $length_flank5 - $length_flank3;
    $ortheus_alignment->{"lengths"} = [$length_flank5, $length_allele, $length_flank3];
    my $ref_allele = substr($ortheus_alignment->{"ref"}->{aligned_sequence}, $length_flank5, $length_allele);
    my $sis_allele = substr($ortheus_alignment->{"sis"}->{aligned_sequence}, $length_flank5, $length_allele);
    my $anc_allele = substr($ortheus_alignment->{"anc"}->{aligned_sequence}, $length_flank5, $length_allele);
    my $old_allele;
    if ($ortheus_alignment->{"old"}) {
        $old_allele = substr($ortheus_alignment->{"old"}->{aligned_sequence}, $length_flank5, $length_allele);
    }
    #    print_ortheus_alignment($ortheus_alignment);
    
    for (my $pos = 0; $pos < length($anc_allele); $pos++) {
        my $anc = substr($anc_allele, $pos, 1);
        my $sis = substr($sis_allele, $pos, 1);
        my $old = "";
        if ($old_allele) {
            $old = substr($old_allele, $pos, 1);
        }
        if ($anc eq $sis and $anc eq $old) {
            $anc = uc($anc)
        } elsif ($anc eq $sis or $anc eq $old) {
            $anc = lc($anc);
        } else {
            $anc = "?";
        }
        substr($anc_allele, $pos, 1, $anc);
    }
    
    $ref_allele =~ s/\-//g;
    $sis_allele =~ s/\-//g;
    $anc_allele =~ s/\-//g;
    $old_allele =~ s/\-//g if ($old_allele);
    
    return($ref_allele, $sis_allele, $anc_allele, $old_allele, $length_flank5, $length_allele, $length_flank3);
}


=head2 print_ortheus_alignment
 
 Arg[1]:        hashref $sorted_alignment
 Returns:       N/A
 Description:   Prints the sorted alignment to the STDOUT. If the alignment has "lengths", this
                method uses those to visually split the alignment in 5' - allele - 3' regions.
                The meothod also prints the "ref", "sis", "anc" and "old" sequences at the end
                if these are defined (these are just references to the corresponding sequences).
                These calls are made by the run_ortheus() method.
                This method is used to print the alignments in verbose mode

=cut

sub print_ortheus_alignment {
    my ($ortheus_alignment) = @_;
    
    foreach my $seq (sort {length($a) <=> length($b) || ($a =~ /^(\d+)/)[0] <=> ($b =~ /^(\d+)/)[0]} grep {!/^ref|sis|anc|old|lengths$/} keys %$ortheus_alignment) {
        if ($ortheus_alignment->{"lengths"}) {
            my $pattern = join("", map {"(.{".$_."})"} @{$ortheus_alignment->{"lengths"}});
            print "SEQ ", join(" | ", ($ortheus_alignment->{$seq}->{aligned_sequence} =~ /$pattern/)), " $seq ", ($ortheus_alignment->{$seq}->{name} or ""), "\n";
        } else {
            print "SEQ ", $ortheus_alignment->{$seq}->{aligned_sequence}, " $seq ", ($ortheus_alignment->{$seq}->{name} or ""), "\n";
        }
    }
    
    foreach my $seq ("ref", "sis", "anc", "old") {
        next if (!$ortheus_alignment->{$seq});
        if ($ortheus_alignment->{"lengths"}) {
            my $pattern = join("", map {"(.{".$_."})"} @{$ortheus_alignment->{"lengths"}});
            print "$seq ", join(" | ", ($ortheus_alignment->{$seq}->{aligned_sequence} =~ /$pattern/)), "\n";
        } else {
            print "$seq ", $ortheus_alignment->{$seq}->{aligned_sequence}, "\n";
        }
    }
}

##############################################################################################
##############################################################################################
##
## SORTED ALIGNMENT METHODS
##
##############################################################################################
##############################################################################################

=head1 Sorted alignment methods

=head2 straighten_sorted_alignment
 
 Arg[1]:        hashref $sorted_alignment
 Returns:       N/A
 Description:   This method reverse-complement the $sorted_alignment if the ref_seq is on the
                reverse strand. I no ref_seq is specified in the $sorted_alignment hashref,
                the first sequence is assumed to be the reference one.
 
=cut

sub straighten_sorted_alignment{
    my ($sorted_alignment) = @_;
    
    my $ref_seq = ($sorted_alignment->{ref_seq} or 0);
    if ($sorted_alignment->{sequences}->[$ref_seq]->{strand} == -1) {
        foreach my $this_seq (@{$sorted_alignment->{sequences}}) {
            if ($this_seq->{aligned_sequence}) {
                $this_seq->{aligned_sequence} = reverse($this_seq->{aligned_sequence});
                $this_seq->{aligned_sequence} =~ tr/ACTGNactgn/TGACNtgacn/;
            }
            if ($this_seq->{original_sequence}) {
                $this_seq->{original_sequence} = reverse($this_seq->{original_sequence});
                $this_seq->{original_sequence} =~ tr/ACTGNactgn/TGACNtgacn/;
            }
            $this_seq->{strand} = -$this_seq->{strand};
        }
    }
}


=head2 get_sub_sorted_alignment
 
 Arg[1]:        hashref $sorted_alignment
 Arg[2]:        int $alignment_offset (first column is 0)
 Arg[3]:        int $alignment_length
 Returns:       hashref $sub_sorted_alignment
 Description:   This method selected a region from a sorted alignment. Starting from
                $alignment_offset and for $alignment_length columns, the resulting
                $sub_sorted_alignment corresponds to that region of the alignment
                only. The coordinates of the sequences and their name are adjusted
                accordingly. The tree labels are not affected.

=cut

sub get_sub_sorted_alignment {
    my ($sorted_alignment, $alignment_offset, $alignment_length) = @_;
    my $sub_sorted_alignment;
    
    $sub_sorted_alignment->{tree} = $sorted_alignment->{tree};
    $sub_sorted_alignment->{positions} = [@{$sorted_alignment->{positions}}];
    
    foreach my $this_seq (@{$sorted_alignment->{sequences}}) {
        my $this_subseq;
        $this_subseq->{species} = $this_seq->{species};
        $this_subseq->{chr} = $this_seq->{chr};
        $this_subseq->{strand} = $this_seq->{strand};
        
        my $flank5 = substr($this_seq->{aligned_sequence}, 0, $alignment_offset);
        $flank5 =~ s/\-//g;
        
        my $aln_seq = substr($this_seq->{aligned_sequence}, $alignment_offset, $alignment_length);
        $aln_seq =~ s/\-//g;
        
        if ($this_subseq->{strand} == 1) {
            $this_subseq->{start} = $this_seq->{start} + length($flank5);
            $this_subseq->{end} = $this_subseq->{start} + length($aln_seq) - 1;
        } elsif ($this_subseq->{strand} == -1) {
            $this_subseq->{end} = $this_seq->{end} - length($flank5);
            $this_subseq->{start} = $this_subseq->{end} - length($aln_seq) + 1;
        } else {
            die "I don't undestand the strand <". $this_subseq->{strand}. ">\n";
        }
        
        $this_subseq->{original_sequence} = $aln_seq; # Still unaligned, but will be aligned by Ortheus.
        
        $this_subseq->{name} = get_name_from_data($this_subseq->{species}, $this_subseq->{chr}, $this_subseq->{start}, $this_subseq->{end}, $this_subseq->{strand});
        
        push(@{$sub_sorted_alignment->{sequences}}, $this_subseq);
    }
    
    return $sub_sorted_alignment;
}


=head2 splice_uninformative_sequences_from_alignment

 Arg[1]:        hashref $sorted_alignment
 Returns:       N/A
 Description:   Look for empty sequences or sequences containing Ns only and splice them
                from the alignment. Note that this method looks at the "original_sequence"
                and not at the "aligned_sequence". This is because this method is used to
                remove uninfromative sequences from the sub-sorted_alignment before the
                sequences are aligned.
                This method is remove the uninformative sequences sequentially until they
                all contain inforamtive nucleotides. The removal of sequences is made via
                the splice_sequence_from_sorted_alignment method.

=cut

sub splice_uninformative_sequences_from_alignment {
    my ($sub_sorted_alignment) = @_;
    my $continue_loop = 0;
#    print $sub_sorted_alignment->{tree}, "\n";
    do {
        $continue_loop = 0;
        for (my $seq_num = 1; $seq_num <= $#{$sub_sorted_alignment->{sequences}}; $seq_num++) {
            my $this_seq = $sub_sorted_alignment->{sequences}->[$seq_num]->{original_sequence};
            if ($this_seq =~ /^[Nn]*$/) {
                $continue_loop = 1;
                splice_sequence_from_sorted_alignment($sub_sorted_alignment, $seq_num);
                last;
            }
        }
#        print $sub_sorted_alignment->{tree}, "\n";
    } while ($continue_loop);
}


=head2 splice_sequence_from_sorted_alignment
 
 Arg[1]:        hashref $sorted_alignment
 Arg[2]:        int $position (starting at 0)
 Returns:       N/A
 Description:   This method splices the sequence in position $position from the input
                $sorted_alignment. It splices the sequence from the "sequences" array,
                adjust the "positions" array and trim the "tree" using the splice_tree
                method.
 
=cut

sub splice_sequence_from_sorted_alignment {
    my ($sorted_alignment, $position) = @_;
    
    # Splice the sequence
    splice(@{$sorted_alignment->{sequences}}, $position, 1);
    
    # Splice the positions array
    for (my $i = 0; $i < @{$sorted_alignment->{positions}}; $i++) {
        if ($sorted_alignment->{positions}->[$i] == $position) {
            splice(@{$sorted_alignment->{positions}}, $i, 1);
            last;
        }
    }
    for (my $i = 0; $i < @{$sorted_alignment->{positions}}; $i++) {
        if ($sorted_alignment->{positions}->[$i] > $position) {
            $sorted_alignment->{positions}->[$i]--;
        }
    }
    
    # Trim the tree
    my $tree = parse_newick($sorted_alignment->{tree});
    $sorted_alignment->{tree} = tree_to_string(splice_tree($tree, $position));
}


##########################################################################################
##########################################################################################
##
## TREE PARSING AND EDITING METHODS
##
##########################################################################################
##########################################################################################

=head1 Tree parsing and editing methods

=cut


=head2 splice_tree
 
 Arg[1]:        $tree
 Arg[2]:        string $leaf_name
 Returns:       $spliced_tree
 Description:   This method splices the leaf $leaf_name from the input $tree.
                It edits the parent node of the removed leaf to avoid returning
                a tree with an internal node that has just one child:

                       o                 o               o
                      / \               / \             / \
                     /   \             /   \           /   \
                    o     o    ==>    o     o   ==>   o     \
                   / \   / \         / \     \       / \     \
                  o   o X   o       o   o     o     o   o     o

                If the leaf to be removed is the only other child of the root
                node, the method will return a tree whose root is the root of
                the remaining tree.

                        o                 o
                       / \                 \                o
                      /   o                 o              / \
                     /   / \     ==>       / \    ==>     o   \
                    /   o   \             o   \          / \   \
                   /   / \   \           / \   \        o   o   o
                  x   o   o   o         o   o   o

                The method works by first calling mark_node() to identify the leaf to be
                removed as well as its parent node. The trim_tree() method is then used to
                edit the tree accordingly.

                This method is usually called from the
                splice_sequence_from_sorted_alignment() method.

=cut

sub splice_tree {
    my ($tree, $leaf_name) = @_;
    
    my $leaf_counter = -1;
    
    mark_node($tree, $leaf_name);
    $tree = trim_tree($tree);
    
    return $tree;
}


=head2 trim_tree

 Arg[1]:        $tree_node
 Returns:       $tree_node
 Description:   This method traverses the tree recursively looking for the node marked
                with a value of 2 (parent node of the leaf to be trimmed, see mark_node()
                and splice_tree()). It assumes that the tree is binary, at least for this
                particular node.
                
                It removes the leaf node to be trimmed and modifies its distance-to-parent
                for its later re-connection to its grand-father node. This is preformed as
                part of the recursivity of the method. The method returns the current node
                except after trimming the parent and leaf node in which case it returns
                the sister node (with its modified distance-to-parent). In each loop, the
                children of the internal nodes are re-connected to the returned value of
                the previous loop, resulting in the expected trimmed tree.

=cut

sub trim_tree {
    ## Marks are:
    ## 1 - leaf node at position $leaf_position
    ## 2 - parent node of leaf node at position $leaf_position
    my ($node) = @_;
    
    if ($node->{mark} and $node->{mark} == 2) {
        if (@{$node->{children}} > 2) {
            die;
        } elsif (@{$node->{children}} == 2) {
            for (my $child_count2 = 0; $child_count2 < 2; $child_count2++) {
                my $this_child = $node->{children}->[$child_count2];
                if (!$this_child->{mark}) {
                    $this_child->{distance_to_parent} += $node->{distance_to_parent};
                    return $this_child;
                }
            }
        } else {
            die;
        }
    }
    
    if ($node->{children}) {
        for (my $i = 0; $i < @{$node->{children}}; $i++) {
            my $this_child = $node->{children}->[$i];
            $node->{children}->[$i] = trim_tree($this_child);
        }
    }
    return $node;
}


=head2 mark_node

 Arg[1]:        $tree_node
 Arg[2]:        string $leaf_name
 Returns:       int $mark (1 if node is leaf, 0 otherwise)
 Description:   This method traverses the tree recursively to find the leaf identified by
                its name. When this happens, the method marks the node with $node->{mark}
                = 1 and returns 1, signalling to the parent node that it is the parent
                node. The parent node is marked with a value of 2.

=cut

sub mark_node {
    ## Marks are:
    ## 1 - leaf node at position $leaf_position
    ## 2 - parent node of leaf node at position $leaf_position
    my ($node, $leaf_name) = @_;
    
    if ($node->{children}) {
        foreach my $this_child (@{$node->{children}}) {
            my $mark = mark_node($this_child, $leaf_name);
            if ($mark == 1) {
                $node->{mark} = 2;
            }
        }
        return 0;
        
    } else {
        if ($node->{label} eq $leaf_name) {
            $node->{mark} = 1;
            return 1;
        } elsif ($node->{label} > $leaf_name) {
            $node->{label}--;
        }
        return 0;
    }
}


=head2 parse_newick

 Arg[1]:        string $newick_tree
 Returns:       hashref $tree_node
 Description:   This method reads the newick tree string recursively and creates a nested
                set of $tree_nodes that represent the tree.

=cut

sub parse_newick {
    my ($newick) = @_;
    my $tree;
    
    
    if ($newick =~ /^\(/) {
        ## This is an internal node. I will go through all the characters, counting opening and closing brakets
        ## to find the commas (,) and last closing bracket corresponding to this internal node.
        
        my $last_comma = 0; # Set marker where the last node separator was
        my $level = 1; # nested level
        
        for (my $i = 1; $i < length($newick); $i++) {
            if (substr($newick, $i, 1) eq "(") {
                ## Additional bracket. Increase level by 1
                $level++;
            } elsif (substr($newick, $i, 1) eq "," and $level == 1) {
                ## Comma at this level. From last_comma to this point, we have a sub-tree:
                push(@{$tree->{children}}, parse_newick(substr($newick, $last_comma + 1, $i-$last_comma - 1)));
                $last_comma = $i;
            } elsif (substr($newick, $i, 1) eq ")") {
                ## Closing bracket. Decrease the level by 1
                $level--;
                if ($level == 0) {
                    ## Closing bracket for this internal node. From last comman to this point we have the last sub-tree:
                    push(@{$tree->{children}}, parse_newick(substr($newick, $last_comma + 1, $i - $last_comma - 1)));
                    ## Assumes no label for internal node, but this could go here if required
                    ## Add distance to parent if exists:
                    if (substr($newick, $i+1) =~ /^\:([\d\.]+)/) {
                        $tree->{distance_to_parent} = $1;
                    } else {
                        $tree->{distance_to_parent} = 0.00001;
                    }
                    return $tree;
                }
            }
        }
        die "I cannot find the closing bracket in the newick string $newick";
    } elsif ($newick =~ /^((\w+)\:([\d\.]+))/) {
        ## Leaf node. Assumes simple label (\w+) and distance (mandatory for this applciation).
        $tree->{label} = $2;
        $tree->{distance_to_parent} = $3;
        return $tree;
    } elsif ($newick =~ /^(\w+)/) {
        ## Leaf node. Assumes simple label (\w+) and distance (mandatory for this applciation).
        $tree->{label} = $1;
        $tree->{distance_to_parent} = 0.0001;
        return $tree;
    } else {
        die "I cannot parse the newick string $newick";
    }
}


=head2 tree_to_string

 Arg[1]:        hashref $tree_node
 Returns:       string $newick_tree
 Description:   This method converts a nested set of tree_nodes into a newick string. It
                calls the node_to_string() method recursively and add a semi-colon (;) at
                the end of the string.

=cut

sub tree_to_string {
    my ($tree) = @_;
    
    my $string = node_to_string($tree). ";";
    
    return $string;
}


=head2 node_to_string

 Arg[1]:        hashref $tree_node
 Returns:       string $newick_node
 Description:   This method writes a Newick representation of this particular node by
                traversing all the children (for internal nodes) or printing the node
                label (for leaves) and prints the distance to the parent node afterwards
                (when available, in principle always but for the root node).

=cut

sub node_to_string {
    my ($node) = @_;
    my $string = "";
    
    if ($node->{children}) {
        $string .= "(";
        $string .= node_to_string($node->{children}->[0]);
        for (my $child_number = 1; $child_number < @{$node->{children}}; $child_number++) {
            $string .= ",";
            $string .= node_to_string($node->{children}->[$child_number]);
        }
        $string .= ")";
    } else {
        $string .= $node->{label};
    }
#    if ($node->{mark}) {
#        $string .= "*" x $node->{mark};
#    }
    if ($node->{distance_to_parent}) {
        $string .= ":".$node->{distance_to_parent};
    }
    return $string;
}




1;