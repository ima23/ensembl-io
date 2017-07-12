=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package Bio::EnsEMBL::IO::Adaptor::VCFAdaptor;
use strict;

use EnsEMBL::Web::Utils::FormatText qw(date_format);
use File::Path qw(make_path);
use Bio::DB::HTS::VCF;

use parent qw(Bio::EnsEMBL::IO::Adaptor::HTSAdaptor);

my $DEBUG = 0;

sub hts_open {
  my $self = shift;

  ## Tabix will want to write the downloaded index file to 
  ## the current working directory. By default this is '/'
  my $time = date_format(time(), '%y-%m-%d'); 
  my $path = $SiteDefs::ENSEMBL_USERDATA_DIR."/temporary/vcf_tabix/$time/";
  make_path($path);
  chdir($path);

  $self->{_cache}->{_htsobj_handle} ||= Bio::DB::HTS::VCF->new(filename => $self->url);
  return $self->{_cache}->{_htsobj_handle};
}

sub htsfile_open {
  my $self = shift;
  
  if (!$self->{_cache}->{_htsfile_handle}) {
    #if (Bio::DB::HTSfile->can('set_udc_defaults')) {
    #  Bio::DB::HTSfile->set_udc_defaults;
    #}
    $self->{_cache}->{_htsfile_handle} = $self->{_cache}->{_htsobj_handle}->{vcf_file}; 
  }
  return $self->{_cache}->{_htsfile_handle};
}

sub htsfile_close {
  my $self = shift;
  
  if ($self->{_cache}->{_htsfile_handle}) {
    $self->{_cache}->{_htsfile_handle}->close();
  }
  return;
}

sub get_header {
  my $self = shift;

  my $hts_obj = $self->hts_open;
  unless ($hts_obj) {
    warn "Failed to open file " . $self->url;
    return undef;
  }
  my $header = $hts_obj->header;
  unless ($header) {
    warn "Failed to get header for file " . $self->url;
    return undef;
  }
  return $header;
}


# UCSC prepend 'chr' on human chr ids. These are in some of the BAM
# files. This method returns a possibly modified chr_id after
# checking whats in the bam file
sub munge_chr_id {
  my ($self, $chr_id) = @_;

  my $header    = $self->get_header;
  return undef unless $header;
  my $seqnames  = $header->get_seqnames || [];

  # Check we get values back for seq region. May need to add 'chr' or 'Chr'
  if (grep { $_ eq $chr_id} @$seqnames) {
    $self->htsfile_close() ;
    return "$chr_id";
  }
  elsif (grep {$_ eq "chr$chr_id"} @$seqnames) {
    $self->htsfile_close() ;
    return "chr$chr_id";
  }
  return undef;
}

sub fetch_variations {
### Get variant info from a VCF file
  my ($self, $chr_id, $start, $end) = @_;

  if ($DEBUG) {
    warn "*** variations for: $chr_id, $start, $end\n";
  }

  my $vcf = $self->hts_open;
  return [] unless $vcf;

  # Maybe need to add 'chr'
  my $seq_id = $self->munge_chr_id($chr_id);
  return [] if !defined($seq_id);

  my @variants = (); # $vcf->get_features_by_location(-seq_id => $seq_id, -start => $start, -end => $end);
  warn ">>> VARIANTS @variants";

  return \@variants;
}


=pod

sub fetch_variations {
  my ($self, $chr, $s, $e) = @_;

  if (!$self->{_cache}->{features} || (ref $self->{_cache}->{features} eq 'ARRAY' && !@{$self->{_cache}->{features}})){
    my @features;
    delete $self->{_cache}->{features};

    ## Eagle fix - tabix will want to write the downloaded index file to 
    ## the current working directory. By default this is '/'
    my $time = date_format(time(), '%y-%m-%d'); 
    my $path = $SiteDefs::ENSEMBL_USERDATA_DIR."/temporary/vcf_tabix/$time/";
    make_path($path);

    foreach my $chr_name ($chr,"chr$chr") { # maybe UCSC-type names?
      my %args = ( 
        region  => "$chr_name:$s-$e",
        file    => $self->url,
        tabix   => $self->hub->species_defs->TABIX,
        tmp_dir => $path,
      );

      my $vcf = Vcf->new(%args);

      ## Eagle fix - this should be called before calling $vcf->next_line to avoid
      ## header line error messages, but require latest vcftools 
      $vcf->parse_header();

      while (my $line=$vcf->next_line()) {
        my $x=$vcf->next_data_hash($line);
        push @features, $x;
      }
      last if(@features);
    }
    $self->{_cache}->{features} = \@features;
  }
  return $self->{_cache}->{features};
}

=cut

1;
