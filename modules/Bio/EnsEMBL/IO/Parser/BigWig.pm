=pod

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

=head1 NAME

Bio::EnsEMBL::IO::Parser::BigWig - A line-based parser devoted to BigWig

=cut

package Bio::EnsEMBL::IO::Parser::BigWig;

use strict;
use warnings;
no warnings 'uninitialized';

use parent qw/Bio::EnsEMBL::IO::BigFileParser Bio::EnsEMBL::IO::Parser::Wig/;
 
=head2 type

    Description : Return case-correct version of format name, for use in method names 
    Returntype  : String

=cut

sub type {
    return 'bigWig'; 
}

=head2 seek

    Description: Fetches the raw data from the requested region and caches it 
    Returntype : Void

=cut

sub seek {
    my ($self, $chr_id, $start, $end) = @_;

    my $fh = $self->open_file;
    warn "Failed to open file ".$self->url unless $fh;
    return unless $fh;

    ## Get the internal chromosome name
    my $seq_id = $self->cache->{'chromosomes'}{$chr_id};
    return unless $seq_id;

    my $intervals = $fh->get_intervals("$seq_id", $start-1, $end);
    my $feature_cache = $self->cache->{'features'};
    foreach my $i (@{$intervals}) {
      push(@{$feature_cache}, [$chr_id, ($i->{start}), ($i->{end}+1), $i->{value}]);
    }
    ## pre-load peek buffer
    $self->next_block();
}

=head2 fetch_summary_data

    Description: fetches data from the requested region, grouped into
                  a set number of bins, and caches it
    Returntype : Void

=cut

sub fetch_summary_data {
    my ($self, $chr_id, $start, $end, $bins) = @_;

    my $fh = $self->open_file;
    warn "Failed to open file ".$self->url unless $fh;
    return unless $fh;

    ## Get the internal chromosome name
    my $seq_id = $self->cache->{'chromosomes'}{$chr_id};
    return unless $seq_id;

    my $list = $fh->get_stats("$seq_id", $start-1, $end, "mean", $bins);
    my $bin_size = floor(($end - $start)/$bins);

    my $feature_cache = [];

    foreach my $value (@$list) {
      next unless defined($_);
      my $line = [$chr_id, $start, $start + $bin_size, $value];
      $start += $bin_size;
      push @$feature_cache, $line;
    }

    $self->cache->{'summary'} = $feature_cache;
}

=head2 fetch_summary_array

    Description: fetches values only from the requested region
    Returntype : ArrayRef

=cut

sub fetch_summary_array {
    my ($self, $chr_id, $start, $end, $bins) = @_;

    my $fh = $self->open_file;
    warn "Failed to open file ".$self->url unless $fh;
    return unless $fh;

    ## Get the internal chromosome name
    my $seq_id = $self->cache->{'chromosomes'}{$chr_id};
    return unless $seq_id;

    ## Get whole chromosome if not defined
    unless ($start && $end) {
      $start = 1;
      $end   = $self->cache->{'chr_sizes'}{$chr_id};
    }

    return $fh->get_stats("$seq_id", $start-1, $end, "mean", $bins);
}

=head2 get_raw_chrom

    Description: Getter for chrom field
    Returntype : String 

=cut

sub get_raw_chrom {
  my $self = shift;
  return $self->{'record'}[0];
}

=head2 get_raw_start

    Description: Getter for start field
    Returntype : Integer 

=cut
sub get_raw_start {
  my $self = shift;
  return $self->{'record'}[1];
}

=head2 get_start

    Description: Getter - wrapper around get_raw_start. Since bigWig features are
                  effectively bedGraph lines, they have semi-open coordinates
    Returntype : Integer 

=cut

sub get_start {
  my $self = shift;
  return $self->get_raw_start() + 1;
}


=head2 get_end

    Description: Getter - wrapper around get_raw_end 
    Returntype : String 

=cut

sub get_end {
  my $self = shift;
  return $self->get_raw_end();
}

=head2 get_raw_score

    Description: Getter for score field
    Returntype : Number (usually floating point) or String (period = no data)

=cut

sub get_raw_score {
  my $self = shift;
  return $self->{'record'}[3];
}


1;
