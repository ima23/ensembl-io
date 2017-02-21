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

BigFileParser - a parser for indexed files such as BigBed and BigWig 

=cut

package Bio::EnsEMBL::IO::BigFileParser;

use strict;
use warnings;

use List::Util qw(max);
use POSIX qw(floor ceil);

use Bio::DB::Big;
use Bio::EnsEMBL::IO::Utils;

use parent qw(Bio::EnsEMBL::IO::Parser);

BEGIN {
  Bio::DB::Big->init();
}

=head2 new

    Constructor
    Argument [1] : URL of file
    Argument [2] : Hash of parameters for configuration, e.g. buffer sizes or 
                   specific functions for handling headers or data
    Returntype   : Bio::EnsEMBL::IO::BigFileParser

=cut

sub new {
    my ($class, $url, $args) = @_;
    
    my $self = {
      url               => $url,
      cache             => {
                            'file_handle' => undef,
                            'features'    => [],
                            },
      current_block     => undef,
      waiting_block     => undef,
      record            => undef,
      strand_conversion => {'+' => '1', '.' => '0', '-' => '-1'},
      %{$args||{}},
    };

    bless $self, $class;
  
    return $self;
}

=head2 open

    Constructor
    Argument [1] : Filepath or GLOB or open filehandle
    Argument [2+]: Hash of parameters for configuration, e.g. buffer sizes or 
                   specific functions for handling headers or data
    Returntype   : Bio::EnsEMBL::IO::BigFileParser

=cut

sub open {
    my ($caller, $url, @other_args) = @_;
    return unless $url;

    ## Trim any whitespace from the URL
    $url =~ s/^\s+//g;
    $url =~ s/\s+$//g;

    my $class = ref($caller) || $caller;
    my $self = $class->new($url, @other_args);

    ## Open and cache the file handle
    my $fh = $self->open_file;
    return unless $fh;
    #warn ">>> OPENED FILE WITH $fh";
 
    ## Cache the chromosome list from the file, mapping Ensembl's non-'chr' names 
    ## to the file's actual chromosome names
    my $chromosomes_hash = $fh->chroms();
    
    my $chromosomes = {};
    my $chr_sizes   = {};
    foreach my $chr_name (keys %{$chromosomes_hash}) {
      my $details = $chromosomes_hash->{$chr_name};
      my $chr = "$chr_name";
      $chr =~ s/^chr//;
      $chromosomes->{$chr} = $details->{name};
      $chromosomes->{$chr_name} = $details->{name};
      $chr_sizes->{$chr} = $details->{length};
    }
    $self->{cache}{chromosomes} = $chromosomes;
    $self->{cache}{chr_sizes}   = $chr_sizes;

    ## Do any additional pre-processing
    $self->init($fh); 

    return $self;
}

sub _query_file {
  my ($self, $chr_id, $callback) = @_;
  my $fh = $self->open_file;
  warn "Failed to open file ".$self->url unless $fh;
  return unless $fh;
  
  ## Get the internal chromosome name
  my $seq_id = $self->cache->{chromosomes}{"$chr_id"};
  if(!$seq_id) {
    warn 'No seq id could be found for '.$chr_id;
    return;
  }
  
  return $callback->($fh, $seq_id);
}

=close

    Description : Dummy method in case a script calls close on a file
    Returntype  : True

=cut

sub close { return 1; }

=head2 init 

    Description: placeholder - may need implementing in child

=cut 

sub init {}


=head2 type

    Description : Placeholder for accessor 
    Returntype  : String

=cut

sub type {
      confess("Method not implemented. This is really important");
}

=head2 url

    Description : Accessor for file url
    Returntype  : String

=cut

sub url {
  my $self = shift;
  return $self->{'url'};
}

=head2 cache

    Description : Accessor for cache
    Returntype  : Hashref

=cut

sub cache {
  my $self = shift;
  return $self->{'cache'};
}


=head2 open_file

    Description: Opens a remote file from URL
    Returntype : Filehandle 

=cut

sub open_file {
  my $self = shift;
  $self->{cache}->{file_handle} ||= Bio::DB::Big->open($self->url);
  return $self->{cache}->{file_handle};
}

=head2 next_block

    Description: Shifts to next block. Note that Big files don't have metadata 
    Returntype : Void

=cut

sub next_block {
    my $self = shift;
    $self->shift_block();
}

=head2 read_block

    Description : Reads a line of text, stores it into next_block, 
                  moving next_block to current_block.
    Returntype   : True/False on existence of a defined current_block after running.

=cut

sub read_block {
    my $self = shift;
    my $features = $self->{'cache'}{'features'};

    if (scalar @$features) {
        $self->{'waiting_block'} = shift @$features || confess ("Error reading cached features: $!");
    } else {
        $self->{'waiting_block'} = undef;
    }
}

=head2 read_record

    Description: Features are cached as an array, so no processing needed 
    Returntype : Void 

=cut


sub read_record {
    my $self = shift;
    $self->{'record'} = $self->{'current_block'};
    #use Data::Dumper; warn '!!! RECORD '.Dumper($self->{'record'});
}


1;
