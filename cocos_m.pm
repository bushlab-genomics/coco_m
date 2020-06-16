=head1 CONTACT

 author: Mariusz Butkiewicz <mariusz.butkiewicz@case.edu>
 modifier & maintainer: BowenJin <bxj139@case.edu>

=cut

=head1 NAME

=head1 SYNOPSIS 

=head1 DESCRIPTION

=cut

package cocos_m;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::Perl;
use List::MoreUtils qw( pairwise );
use List::MoreUtils qw( each_array );
use feature qw(say);

my $termination_code = 0;
my $stop_codon_pos_5prime = 0; 

sub version {
    return '80';
}

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {
        ALT_AA_SEQ => "alterated amino acid sequence (COCOS Plugin)",
        TERM_TYPE => "termination type of sequence(COCOS Plugin)",     
	IMPACT => "impact of variant (COCOS Plugin)",
	PTC_POSITION => "the new premature codon position (COCOS Plugin)"
    };
}


#sub new {
#    my $class = shift;
#    my $self = $class->SUPER::new(@_);
#    
#    # configure
#    my @params = @{$self->params};
#    $self->{fasta}     = $params[0] || "./cocos.fasta";
#    return $self;
#}

sub run {
    my ($self, $transcript_variation_allele, $vep_line) = @_;

    my $trans = $transcript_variation_allele->transcript;
    my $variant_cdna_pos = $vep_line->{'cDNA_position'};
    my $variant_cds_pos = $vep_line->{'CDS_position'};
    my $transcript_id = $vep_line->{'Feature'};
    my $uploaded_variation = $vep_line->{'Uploaded_variation'};
    my $consequence = $vep_line->{'Consequence'};
    $uploaded_variation =~ s/\//_/g;

    return {} if $consequence !~ /stop_lost/ && $consequence !~ /frame/ || $consequence =~ /NMD/;
	
    if($variant_cdna_pos !~ /^-$/ && $variant_cds_pos !~ /^-$/) {

        $variant_cdna_pos = (split /-/, $variant_cdna_pos)[0]; 
	
        my $result = process_transcript($transcript_variation_allele, $variant_cdna_pos);
	my $result_aa = translate_seq_string($result);

	my $stop_csq = '';
	my $stop_impact = '';
	if ($termination_code == 1) { 
		$stop_csq = "PTC, NMD";
		$stop_impact = "mild";
	} elsif ($termination_code == 2) {
		$stop_csq = "NO PTC, NMD";
		$stop_impact = "mild";
	} elsif ($termination_code == 3) { 
		$stop_csq = "PTC,canonical,escaping NMD";
		$stop_impact = "high severity";
	} elsif ($termination_code == 4) { 
		$stop_csq = "PTC,noncanonical,escaping NMD";
		$stop_impact = "mild severity";
	} else { ; }

        return {
            ALT_AA_SEQ => $result_aa,
	    TERM_TYPE => $stop_csq,
	    IMPACT => $stop_impact,
	    PTC_POSITION => $stop_codon_pos_5prime
        };
    }

    return {};
}

# translate a DNA sequence string into a amino acid sequence string
# if sequence is not empty or contains 'X' characters (unknown)
sub translate_seq_string {
    my ($custom_seq) = @_;

    if($custom_seq =~ /X/) {
	print "sequence can not be translated! contains X characters!\n";
	print "$custom_seq\n";
        return "";
    }

    if($custom_seq ne '') {
        my $seq_object = new_sequence($custom_seq,'','');
        return $seq_object->translate->seq();
    }
    return '';
}

#get coding seq excluding 5'UTR but including 3'UTR
sub transcript_seq_spliced_without_5prime_utr {
   my ($tr) = @_;
   my $ori_sequence = $tr->spliced_seq("soft_mask");
   my $sequence = substr $ori_sequence, $tr->cdna_coding_start()-1, length $tr->spliced_seq("soft_mask");
   my $five_prime_utr_len = (length $ori_sequence ) - (length $sequence) ;
   return $sequence, $five_prime_utr_len;
}

#get transcript coding region
sub transcript_coding_seq {
   my ($tr) = @_;
   my $sequence = $tr->spliced_seq("soft_mask");
   $sequence =~ s/[a-z]//g; 
   return $sequence;
}

# slice sequence and insert variant of interest 
sub insert_variant_into_exon {
    my ($sequence, $variant, $tva) = @_;

    my $tr = $tva->transcript;
    my $exon_start = $tr->cdna_coding_start();
    my $exon_end = length $tr->spliced_seq("soft_mask");
    my $var_start = $tva->transcript_variation->cdna_start;
    my $var_end = $tva->transcript_variation->cdna_end;

    my $three_prime_utr_overhead = length($sequence) - ($exon_end - $exon_start + 1);
    return "" if not defined($var_end);

    $variant =~ s/-//;
    my $before_variant = substr($sequence,0,$var_start - $exon_start);
    my $after_variant = substr($sequence, $var_end - $exon_start + 1, $exon_end - $var_end + $three_prime_utr_overhead);

    return ($before_variant . $variant . $after_variant, ($var_start - $var_end - 1) );
}


sub get_last_ejunct {
        my ($trans, $var_cdna_pos, $variant_len) = @_;
        my $ori_ejunct_pos = ("");
        my $var_ejunct_pos = ("");
	my @ori_exon_start_pos = "";
	my @ori_exon_end_pos = "";
	my @var_exon_start_pos = "";
	my @var_exon_end_pos = "";

        my $trans_seq = $trans->spliced_seq;
        my @coding_coord = ($trans->cdna_coding_start(),$trans->cdna_coding_end());


        if(defined $coding_coord[0] && defined $coding_coord[1]){
                my @all_exon = @{$trans->get_all_Exons()};
                @ori_exon_end_pos   = map{$_->cdna_end($trans)}(@all_exon);
                @ori_exon_start_pos = map{$_->cdna_start($trans)}(@all_exon);
                @var_exon_end_pos   = @ori_exon_end_pos;
		@var_exon_start_pos = @ori_exon_start_pos;

                my $idx = 0;
                while($idx < (@ori_exon_end_pos-1) && $ori_exon_end_pos[$idx] < $coding_coord[1]){
                        $idx+=1;
                }

		$ori_ejunct_pos = $ori_exon_end_pos[$idx-1];

		my $var_exon_idx = "";
		my $flag = 0;
		for (my $i=0; $i<=$idx; $i++){
			if ($ori_exon_start_pos[$i] < $var_cdna_pos && $ori_exon_end_pos[$i] > $var_cdna_pos){
				$flag = 1;
				$var_exon_idx = $i;
				$var_exon_end_pos[$i] = $var_exon_end_pos[$i] + $variant_len;
				next;}

			if ($flag){
				$var_exon_start_pos[$i] += $variant_len;
				$var_exon_end_pos[$i]   += $variant_len;}
		}
		$var_ejunct_pos = $var_exon_end_pos[$idx-1]; }

	return ($var_ejunct_pos, @var_exon_end_pos);
}

#split seq in tri-grams denoting a list of codons
sub codon_ize_sequence {
    my ($sequence) = @_;
    $sequence =~ s/(\w{3})/$1 /g;
    return split(/ /,$sequence);
}

#determine whether condon is a stop-codon
sub is_stop_codon {
    my ($codon) = @_;
    $codon = uc $codon;
    return ($codon =~ /TAG/) || ($codon =~ /TAA/) || ($codon =~ /TGA/);
}

sub compare_codon_lists {
    my @codons_w_variant = @{$_[0]};
    my @codons_wo_variant = @{$_[1]};

    my $truncate_sequence = "";
    # stop codon starts
    my $stop_codon_pos = 0;

    my $i = 0;
    for ($i=0; $i<=$#codons_w_variant; $i++ ){
	$truncate_sequence .= $codons_w_variant[$i];
	if (is_stop_codon( $codons_w_variant[$i] )){
		last;}
    }

    if ($i < $#codons_w_variant){
	$stop_codon_pos = $i * 3;
	$termination_code = 1;
    } else {
	$termination_code = 2;
    }
   
    # final viable captured sequence
    return ($truncate_sequence, $stop_codon_pos);
}


sub process_transcript {
    my ($transcript_variation_allele, $variant_cdna_pos) = @_;
    my $transcript = $transcript_variation_allele->transcript();
    my ($seq, $five_prime_UTR_len) = transcript_seq_spliced_without_5prime_utr($transcript);
    my $ref_variant_seq = (split /\//, $transcript_variation_allele->allele_string)[0];
    my $variant = $transcript_variation_allele->variation_feature_seq;

    my ($condon_seq_with_variant, $variant_len) = insert_variant_into_exon
    (
        $seq,
        $variant,
        $transcript_variation_allele
    );

    my($var_ejunct, @var_exon_end_pos) = get_last_ejunct($transcript, $variant_cdna_pos, $variant_len);

    my $condon_seq_wo_variant = $seq;

    my @codons_w_var = codon_ize_sequence($condon_seq_with_variant);
    my @codons_wo_var = codon_ize_sequence($condon_seq_wo_variant);
    my ($truncate_sequence, $stop_codon_pos) = compare_codon_lists(\@codons_w_var, \@codons_wo_var);

    # the distance between EJ and PTC
    $stop_codon_pos_5prime = $stop_codon_pos + $five_prime_UTR_len + 1; 
    my $diff = $var_ejunct - $stop_codon_pos - $five_prime_UTR_len;

    my $PTC_exon_len = 0;
    # the length of the exon where PTC is
    for (my $j=0; $j<=$#var_exon_end_pos; $j++){
      if ($stop_codon_pos < $var_exon_end_pos[$j] - $five_prime_UTR_len){
        $PTC_exon_len = $var_exon_end_pos[$j] - $var_exon_end_pos[$j-1];
      	last;}
    }

    ### canonical rules: 55 nt rule && last exon rule
    if ($stop_codon_pos > 0 && $termination_code == 1 &&  $diff < 55){
	$termination_code = 3;
	return $truncate_sequence;
    ### uncanonical rules
    } elsif ($PTC_exon_len > 407) {
	$termination_code = 4;
	return $truncate_sequence;
    } elsif ($stop_codon_pos < 150) {
	$termination_code = 4;
	return $truncate_sequence;
    } else {
	return "";
    }
}

1;

