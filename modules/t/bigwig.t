# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;

use Test::More;
use Bio::EnsEMBL::IO::Parser::BigWig;

######################################################
## Test 1
######################################################
diag 'Testing variable step bigwig';
my $parser = Bio::EnsEMBL::IO::Parser::BigWig->open("modules/t/input/data-variableStep.bw");
ok($parser->seek(1, 1, 15), 'Can query for the right region (chr 1, start 1, end 15)');

ok($parser->next, 'Can iterate');
is($parser->get_seqname, '1', 'Seq name');
is($parser->get_start, 1, 'Start');
is($parser->get_end, 2, 'End');
is($parser->get_score, 1, 'Score');

for (my $i = 1; $i <= 4; $i++) {
  ok($parser->next, 'Iterate block '.$i);
  is($parser->get_seqname, '1', 'Seq name block '.$i);
  my $start = ($i * 2);
  is($parser->get_start, $start, 'Start '.$i);
  is($parser->get_end, $start + 1, 'End '.$i);
  is($parser->get_score, $i + 1, 'Score '.$i);
}

ok(!$parser->next, 'No more data left');

diag 'Testing summary code';
{  
  my $summary_array = $parser->fetch_summary_array('chr1', 1, 15, 5);
  is(scalar(@{$summary_array}), 5, 'Expect 5 bins returned');
  is_deeply($summary_array,[1.5,3.5,5.0,undef,undef], 'Checking summary array values');


  my $summary_data = $parser->fetch_summary_data('chr1', 1, 15, 5);
  is(scalar(@{$summary_data}), 5, 'Expect 5 bins returned');
  is_deeply($summary_data,[
    ['chr1', 1, 3, 1.5],
    ['chr1', 4, 6, 3.5],
    ['chr1', 7, 9, 5.0],
    ['chr1', 10, 12, undef],
    ['chr1', 13, 15, undef],
  ], 'Checking all bin summaries returned as expected') or diag explain $summary_data;
}

$parser->close();

######################################################
## Test 2
######################################################
diag 'Testing fixed step bigwig';
$parser = Bio::EnsEMBL::IO::Parser::BigWig->open('modules/t/input/data-fixedStep.bw');
ok($parser->seek(1, 1, 15), 'Can query for new file across right region (chr 1, start 1, end 15)');

for (my $i = 1; $i <= 10; $i ++) {
  ok($parser->next, 'Can interate');
  is($parser->get_seqname, '1', 'Seq name '.$i);
  is($parser->get_start, $i, 'Start '.$i);
  is($parser->get_end, $i + 1, 'End '.$i);
  is($parser->get_score, ($i-1), 'Score '.$i);
}

ok(!$parser->next, 'No more data left');
$parser->close();

diag 'Testing summary code';
{  
  my $summary_array = $parser->fetch_summary_array('chr1', 1, 15, 5);
  is(scalar(@{$summary_array}), 5, 'Expect 5 bins returned');
  is_deeply($summary_array,[1.0,4.0,7.0,9.0,undef], 'Checking summary array values');

  my $summary_data = $parser->fetch_summary_data('chr1', 1, 15, 5);
  is(scalar(@{$summary_data}), 5, 'Expect 5 bins returned');
  is_deeply($summary_data,[
    ['chr1', 1, 3, 1.0],
    ['chr1', 4, 6, 4.0],
    ['chr1', 7, 9, 7.0],
    ['chr1', 10, 12, 9.0],
    ['chr1', 13, 15, undef],
  ], 'Checking all bin summaries returned as expected') or diag explain $summary_data;
}

done_testing;
