#!/usr/bin/env perl
my $VERSION = '2.1';
my $lastmodif = '2021-3-5';

use strict;
use warnings;
use diagnostics;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Getopt::Long;
use Data::Dumper;

my $idTag = "id";
my $sourceTag = basename($0) . '-' . $VERSION;
my $opt_seqId;
my $opt_digit = 6;
my $opt_length = 1;
my $opt_note;
my $opt_outfile;
my $opt_RMclariTE;
my $opt_featurePrefix = 'genome';
my $opt_cds;
&GetOptions("h|help"   => \&help,
            "source=s" => \$sourceTag,
            "seqid=s"  => \$opt_seqId,
            "RMclariTE"    => \$opt_RMclariTE,
			"featurePrefix=s" => \$opt_featurePrefix,
            "cds"      => \$opt_cds,
            "o"        => \$opt_outfile,
            "d=i"      => \$opt_digit,
            "note"     => \$opt_note, #convert all tags to a single Note="tag:xxx; tag:xxx; tag:xxx"
            "l=i"      => \$opt_length,
            "id=s"     => \$idTag);

@ARGV or &help;
&main();

sub main {
	my $self = {};
	bless $self;
	$self->setOptions();
	foreach (@ARGV){
		$self->{inputFile} = $_;
		-e $self->{inputFile} or die("file " . $self->{inputFile} . " does not exist");
		$self->createGffIo();
		$self->getFeatures();
	}
	exit(0);
}

sub setOptions {
	my $self = shift;
	$self->{idTag} = $idTag;
	$self->{sourceTag} = $sourceTag;
	$self->{opt_seqId} = $opt_seqId;
	$self->{opt_outfile} = $opt_outfile;
	$self->{opt_digit} = $opt_digit;
	$self->{opt_RMclariTE} = $opt_RMclariTE;
	$self->{opt_featurePrefix} = $opt_featurePrefix;
	$self->{opt_cds} = $opt_cds;
	$self->{opt_note} = $opt_note;
	$self->{opt_length} = $opt_length;
	return 1;
}

sub createGffIo {
	my $self = shift;
	if($self->{opt_outfile}){
		$self->{outFile} = $self->{inputFile} . ".gff";
		$self->{gffIo} = Bio::Tools::GFF->new(-file => ">".$self->{outFile},
		                                      -gff_version => "3");
	}
	else{
		$self->{gffIo} = Bio::Tools::GFF->new(-fh => \*STDOUT,
		                                      -gff_version => "3");
	}
	return 1;
}

sub getFeatures {
	my $self = shift;
	$self->{seqIo} = Bio::SeqIO->new(-format => "EMBL", -file => $self->{inputFile});
	my $kSeq=0;
	while($self->{seqObj} = $self->{seqIo}->next_seq()){
		$kSeq++;
		$self->{seqId} = $self->{seqObj}->display_id() ;
		if(lc($self->{seqId}) eq "unknown"){
			$self->{seqId} = basename($self->{inputFile}) . "." . $kSeq;
		}

		$self->{opt_RMclariTE} and $self->RMclariTE_clarite_seqid_builder();

		### the "region" feature
		$self->{feature}->{$kSeq}->{region} = Bio::SeqFeature::Generic->new(-seq_id      => $self->{seqId},
		                                                                    -source_tag  => $self->{sourceTag},
		                                                                    -primary_tag => 'region',
		                                                                    -start       => 1,
		                                                                    -end         => $self->{seqObj}->length(),
		                                                                    -strand      => 1);
		$self->{gffIo}->write_feature( $self->{feature}->{$kSeq}->{region} );
		
		foreach my $feature ($self->{seqObj}->all_SeqFeatures) {
			### 1st: create the parent feature ###
			$self->correct_locs($feature->location);
			
			#skip small features
			if($self->{feat}->{location}->{end} - $self->{feat}->{location}->{start} <= $self->{opt_length}){
				next;
			}
			$self->{feat}->{obj} = Bio::SeqFeature::Generic -> new ( -seq_id      => $self->{seqId},
			                                                         -source_tag  => $self->{sourceTag},
			                                                         -primary_tag => $feature->primary_tag,
			                                                         -start       => $self->{feat}->{location}->{start},
			                                                         -end         => $self->{feat}->{location}->{end},
			                                                         -strand      => $self->{feat}->{location}->{strand} ) ;

			if($self->{opt_RMclariTE}){
				$self->RMclariTE_clarite_gff_builder($feature);
				$self->{gffIo}->write_feature($self->{feat}->{obj});
				$self->create_subfeatures_match_part();
			}
			elsif($self->{opt_cds}){
				$feature->primary_tag eq "CDS" or next;
				$self->generic_gff_builder($feature);
				$self->create_features_gene_mrna_exon();
			}
			else{
				$self->generic_gff_builder($feature);
				$self->{gffIo}->write_feature($self->{feat}->{obj});
				$self->create_subfeatures_match_part();
			}
		}
	}
	return 1;
}

sub RMclariTE_clarite_seqid_builder {
	my $self = shift;
	if (basename($self->{inputFile})=~/^((\w+)#\d+#\d+)\.fa/) {
		($self->{seqId}, $self->{chr}) = ($1, $2);
	}
	elsif (basename($self->{inputFile})=~/^((\w+):\d+\-\d+)\.fa/) {
		($self->{seqId}, $self->{chr}) = ($1, $2);
	}
	return 1;
}

sub RMclariTE_clarite_gff_builder {
	my $self = shift;
	my $feature = shift;
	my $pritag = $self->{feat}->{obj}->primary_tag;
	$self->{k_pritag}->{$pritag} or $self->{k_pritag}->{$pritag} = 0; #initialize pritag counter
	$self->{k_pritag}->{$pritag}++;
	my $ftid = $self->{opt_featurePrefix} . $self->{chr} . '_' . sprintf $pritag."_%.".$self->{opt_digit}."d", $self->{k_pritag}->{$pritag};
	$self->{feat}->{obj}->add_tag_value("ID", $ftid);
	$self->{feat}->{id} = $ftid;

	### Parse CLARITE feature tags
	my $parsed_tags = $self->parse_tagvalues_clarite($feature);
	foreach my $tag (sort keys %{$parsed_tags}) {
		foreach (sort {$a<=>$b} keys %{$parsed_tags->{$tag}}) {
			$self->{feat}->{obj}->add_tag_value($tag, $parsed_tags->{$tag}->{$_}); 
		}
	}
	return 1;
}

sub parse_tagvalues_clarite {
	my $self = shift;
	my $feature = shift;
	my @tags = $feature->all_tags;
	my $hash_tag;
	foreach my $tag (@tags) {
		lc$tag eq 'id' and next ;
		lc$tag eq 'parent' and next;
		lc$tag eq 'range' and next;
		my $gfftag = uc(substr($tag,0,1)) . substr($tag,1);
		if($tag eq 'compo'){
			($hash_tag->{Family}->{1}, $hash_tag->{Compo}->{1}) = parse_tag_compo($feature); 
		}
		elsif ($tag eq 'copie'){
			 $hash_tag->{Matching_repeat}->{1} = join(" ", $feature ->each_tag_value($tag));
		}
		else{
			my $k=1;
			foreach my $tagvalue ($feature->each_tag_value($tag)) {
				$tagvalue eq "_no_value" and $tagvalue = 'true' ; # This "_no_value" is return by Bio::Seq when you have for example FT     /pseudo
				$tagvalue =~s/^\s+|\s+$//g ;
				$hash_tag->{$gfftag}->{$k++} = $tagvalue;
			}
		}
	}
	# push all tags in 1 Note 
	if ($self->{opt_note}) {
		my @note = ();
		foreach my $tag (sort keys %{$hash_tag}) {
			foreach (sort {$a<=>$b} keys %{$hash_tag->{$tag}}) {
				push @note, $tag.':'.$hash_tag->{$tag}->{$_};
			}
		}
		undef $hash_tag; 
		$hash_tag->{Note}->{1} = join(" ", @note);
	}
	return $hash_tag;
}

sub parse_tag_compo {
	my $feature = shift;
	my $compo_value = ($feature ->each_tag_value("compo"))[0];
	$compo_value =~s/\s+$//;
	my @compo = split/\s+/, $compo_value;
	my @compo_ordered;
	my $family = 'unknown';
	for (my $i=0 ; $i<$#compo ; $i+=2){
		next if($compo[$i] eq 'no_match') ;
		push @compo_ordered, [$compo[$i], $compo[$i+1]];
	}
	if(@compo_ordered){
		@compo_ordered = sort {$b->[1] <=> $a->[1]} @compo_ordered;
		$family = $compo_ordered[0]->[0];
	}
	return($family, $compo_value);
}

sub generic_gff_builder {
	my $self = shift;
	my $feature = shift;
	my $pritag = $self->{feat}->{obj}->primary_tag;
	$self->{k_pritag}->{$pritag} or $self->{k_pritag}->{$pritag} = 0; #initialize pritag counter
	$self->{k_pritag}->{$pritag}++;
	
	my $ftid;
	if( $feature->has_tag($self->{idTag}) ){
		$ftid = ($feature->each_tag_value($self->{idTag}))[0];
	}
	else{
		$ftid = $pritag . '_' . $self->{k_pritag}->{$pritag};
	}

	$self->{feat}->{obj}->add_tag_value("ID", $ftid);
	$self->{feat}->{id} = $ftid;

	### Parse feature tags
	my $parsed_tags = $self->parse_tagvalues_generic($feature);
	foreach my $tag (sort keys %{$parsed_tags}) {
		foreach (sort {$a<=>$b} keys %{$parsed_tags->{$tag}}) {
			$self->{feat}->{obj}->add_tag_value($tag, $parsed_tags->{$tag}->{$_}); 
		}
	}
	return 1;
}

sub parse_tagvalues_generic {
	my $self = shift;
	my $feature = shift;
	my @tags = $feature->all_tags;
	my $hash_tag;
	foreach my $tag (@tags) {
		$tag eq $self->{idTag} and next;
		$tag eq 'translation' and next;
		my $gfftag = uc(substr($tag,0,1)) . substr($tag,1);
		my $k=1;
		foreach my $tagvalue ($feature->each_tag_value($tag)) {
			$tagvalue eq "_no_value" and $tagvalue = 'true' ; # This "_no_value" is return by Bio::Seq when you have for example FT     /pseudo
			$tagvalue =~s/^\s+|\s+$//g ;
			$tagvalue =~s/%/percent/g ;
			$tagvalue =~s/, / /g ;
			$tagvalue =~s/,/ /g ;
			$hash_tag->{$gfftag}->{$k++} = $tagvalue;
		}
	}
	return $hash_tag;
}

sub create_subfeatures_match_part {
	my $self = shift;
	my $location = $self->{feat}->{location}->{all_locs};
	my $ftid = $self->{feat}->{id}; 
	my $subk = 1;
	foreach (@{$location}) {
		my $feat = Bio::SeqFeature::Generic ->new( -seq_id      => $self->{seqId},
		                                           -source_tag  => $self->{sourceTag},
		                                           -primary_tag => "match_part",
		                                           -start       => $_->[0],
		                                           -end         => $_->[1],
		                                           -strand      => $_->[2] );
		$feat ->add_tag_value ( "Parent", $ftid );
		$feat ->add_tag_value ( "ID", $ftid . "_" . $subk++);
		$self->{gffIo}->write_feature($feat);
	} 
	return 1;
}

sub create_features_gene_mrna_exon {
	my $self = shift;
	my $location = $self->{feat}->{location}->{all_locs};
	
	### gene / mRNA
	foreach my $type ('gene', 'mRNA'){
		my $feat = Bio::SeqFeature::Generic ->new( -seq_id      => $self->{seqId},
	        	                                   -source_tag  => $self->{sourceTag},
	        	                                   -primary_tag => $type,
	        	                                   -start       => $self->{feat}->{location}->{start},
	        	                                   -end         => $self->{feat}->{location}->{end},
	        	                                   -strand      => $self->{feat}->{location}->{strand} );
		$feat ->add_tag_value('ID', $self->{feat}->{id} . '_' . $type);
		if($type eq 'mRNA'){
			$feat ->add_tag_value("Parent", $self->{feat}->{id} . '_gene');
			$self->{feat}->{obj}->add_tag_value("Parent", $self->{feat}->{id} . '_mRNA');
		}
		$self->{gffIo}->write_feature($feat);
	}

	### CDS
	$self->{gffIo}->write_feature($self->{feat}->{obj});

	### exons
	my $subk = 1;
	foreach (@{$location}) {
		my $feat = Bio::SeqFeature::Generic ->new( -seq_id      => $self->{seqId},
		                                           -source_tag  => $self->{sourceTag},
		                                           -primary_tag => 'exon',
		                                           -start       => $_->[0],
		                                           -end         => $_->[1],
		                                           -strand      => $_->[2] );
		$feat ->add_tag_value("Parent", $self->{feat}->{id} . '_mRNA');
		$feat ->add_tag_value("ID", $self->{feat}->{id} . "_exon_" . $subk++);
		$self->{gffIo}->write_feature($feat);
	} 
	return 1;
}

sub correct_locs {
	my $self = shift;
	my $locationobj = shift;
	my @locs = sort {$a->start <=> $b->end} $locationobj->each_Location();
	$self->{feat}->{location}->{start} = $locs[0]->start;
	$self->{feat}->{location}->{end} = $locs[$#locs]->end;
	$self->{feat}->{location}->{strand} = $locs[0]->strand;

	my $loc;
	foreach (@locs){
		push @{$loc}, [$_->start, $_->end, $_->strand];
	}
	$self->{feat}->{location}->{all_locs} = $loc;
	return 1;
}

sub help {
	my $prog = basename($0);
	print STDERR <<EOF;
#### $prog $VERSION ####
#
# AUTHOR:     Frederic Choulet
# VERSION:    $VERSION
# LAST MODIF: $lastmodif
# This script is used to print GFF3 files from EMBL files

USAGE:
      $prog  [OPTIONS]  embl_file1  embl_file2  ...

         ### OPTIONS ###
         -id <string>          ID tag in EMBL file [default: /id]
         -seqid <string>       seq_id tag in GFF file [default: the EMBL ID]
                               if EMBL ID = unknown -> seqid = file name
         -source <string>      Source tag in GFF file [default: EMBL2GFF]
         -o                    Create output file [default: STDOUT]
         -d <integer>          Length of the digit in the ID of features [default: 6]
         -note                 Convert all tags/values to a single tag Note="tag:xxx; tag:xxx; tag:xxx"
         -cds                  Build features gene/mRNA/exon according to CDS (considers only CDSs)
         -RMclariTE            Build GFF specifically for RepeatMasker CLARITE EMBL files
		 -featurePrefix        Give custom (genome specific) prefix for CLARITE features, only if RMclariTE option enabled (default: genome)
         -l <interger>         Length threshold to skip small features [default: 1 bp]
         -h:                   Print this help

EOF
	exit(1);
}


