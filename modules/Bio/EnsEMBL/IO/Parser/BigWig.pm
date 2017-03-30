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

=head1 DESCRIPTION

This interacts with the L<Bio::DB::Big> library and handles the translation of 
coordinates from Ensembl style C<1-start, fully-closed> (where indexing starts at 1) to
the UCSC Big coordiante system C<0-start, half-open> (where indexing starts at 0) and 
back again.

=cut

package Bio::EnsEMBL::IO::Parser::BigWig;

use strict;
use warnings;
no warnings 'uninitialized';
use POSIX qw/floor/;

use parent qw/Bio::EnsEMBL::IO::BigFileParser Bio::EnsEMBL::IO::Parser::Wig/;
 
=head2 type

    Description : Return case-correct version of format name
    Returntype  : String

=cut

sub type {
    return 'bigWig'; 
}

=head2 seek

    Description: Fetches the raw data from the requested region and caches it. Use
                 the accessor methods such as <get_start>, <get_end> and <get_score>
                 to extract the information from an iteration
    Returntype : Boolean; true if there is a score to retrieve

=cut

sub seek {
  my ($self, $chr_id, $start, $end) = @_;
  return $self->_query_file($chr_id, sub {
    my ($fh, $seq_id) = @_;
    my $intervals = $fh->get_intervals("$seq_id", $start-1, $end);
    my $feature_cache = [];
    foreach my $i (@{$intervals}) {
      push(@{$feature_cache}, [$chr_id, ($i->{start}), ($i->{end}+1), $i->{value}]);
    }
    $self->cache->{'features'} = $feature_cache;
    ## pre-load peek buffer
    return $self->next_block();
  });
}

=head2 fetch_summary_data

    Description: fetches data from the requested region, grouped into
                  a set number of bins with start and ends. Bins without a value 
                  will still be returned but as an undefined value
    Returntype : ArrayRef[ArrayRef[chromosome, start, end, float_value]]

=cut

sub fetch_summary_data {
  my ($self, $chr_id, $start, $end, $bins) = @_;
  
  my $values = $self->fetch_summary_array($chr_id, $start, $end, $bins);
  # length of requested region / bins floored
  my $bin_size = floor(($end - $start+1)/$bins);
  
  my $features = [];
  foreach my $value (@$values) {
    #If start = 1, end = 30, bins = 3, binsize = 10 = (lengthofregion/bins) = ((30-1+1)/3)
    #iter 1: start = 1.             end = (1+10-1) = 10
    #iter 2: start = (1+bin) =  11. end = (11+10-1) = 20
    #iter 3: start = (11+bin) = 21. end = (21+10-1) = 30
    my $line = [$chr_id, $start, ($start + $bin_size - 1), $value];
    push(@{$features}, $line);
    $start += $bin_size;
  }
  
  return $features;
}

=head2 fetch_summary_array

    Description: fetches values only from the requested region. Calculated by bins
    Returntype : ArrayRef[float_value]

=cut

sub fetch_summary_array {
  my ($self, $chr_id, $start, $end, $bins) = @_;
  return $self->_query_file($chr_id, sub {
    my ($fh, $seq_id) = @_;    
    $start = 1 unless $start;
    $end = $self->cache->{chr_sizes}->{$chr_id} unless $end;
    return $fh->get_stats("$seq_id", $start-1, $end, $bins, "mean");
  });
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

    Description: Getter for start field (0 based)
    Returntype : Integer 

=cut
sub get_raw_start {
  my $self = shift;
  return $self->{'record'}[1];
}

=head2 get_start

    Description: Getter - wrapper around get_raw_start converting into Ensembl style coordinates
    Returntype : Integer 

=cut

sub get_start {
  my $self = shift;
  return $self->get_raw_start() + 1;
}


=head2 get_end

    Description: Getter - returns back the end
    Returntype : String 

=cut

sub get_end {
  my $self = shift;
  return $self->{record}->[2];
}

=head2 get_raw_score

    Description: Getter for score field
    Returntype : Float

=cut

sub get_raw_score {
  my $self = shift;
  return $self->{record}->[3];
}

=head2 get_score

    Description: Getter for score field
    Returntype : Float

=cut

sub get_score {
  my $self = shift;
  return $self->{record}->[3];  
}

1;
