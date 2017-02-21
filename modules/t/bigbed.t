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
use Bio::EnsEMBL::IO::Parser::BigBed;

######################################################
## Test 1
######################################################
my $parser = Bio::EnsEMBL::IO::Parser::BigBed->open('modules/t/input/data.bb');
ok($parser->seek(1, 1, 10), "Seek to location");

ok($parser->next(), 'Fetch next set of rows');
is($parser->get_chrom, '1', 'Check chromosome');
is($parser->get_start, 3, 'Check start');
is($parser->get_end, 6, 'Check end');
is($parser->get_strand, 0, 'Check strand');
is($parser->get_name, 'Mo', 'Check name');
is($parser->get_score, 1000, 'Check score');

ok($parser->next(), 'Fetch next set of rows');
is($parser->get_chrom , '1', 'Check chromosome');
is($parser->get_start , 4, 'Check start');
is($parser->get_end , 8, 'Check end');
is($parser->get_strand, 1, 'Check strand');
is($parser->get_name , 'Larry', 'Check name');
is($parser->get_score , 1000, 'Check score');

ok($parser->seek(2, 1, 10), "Seek to new location");
ok($parser->next(), 'Fetch next set of rows');
is($parser->get_chrom , '2', 'Check chromosome');
is($parser->get_start , 2, 'Check start');
is($parser->get_end , 7, 'Check end');
is($parser->get_strand , -1, 'Check strand');
is($parser->get_name , 'Curly', 'Check name');
is($parser->get_score , 1000, 'Check score');

ok(!$parser->next, 'No more rows');

$parser->close();

done_testing;
