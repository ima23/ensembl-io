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

$parser->close();

######################################################
## Test 2
######################################################
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

done_testing;
