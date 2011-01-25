#!/usr/bin/ruby

require 'csv'
require 'set'
require 'ftools'
require 'optparse'

# get the files from the same directory as us
$: << File.dirname(__FILE__) 

require 'read_file.rb'
require 'protein.rb'
require 'xme_filter.rb'
require 'rename.rb'

$options = {
	:verbose => true,
	:fatal => false,
	:output => "output.csv",
	:output2 => "output2.csv"
}
OptionParser.new do |opts|
	opts.banner = "Usage: #{$0} [options] NP4Row10feature.csv ..."

	opts.on("-o", "--output=FILE", "output to FILE") do |file|
		$options[:output] = file
	end
	opts.on("-f", "--fatal", "Stop on first error") do |v|
		$options[:fatal] = v
	end
	opts.on("-q", "--quiet", "Run quietly") do |v|
		$options[:verbose] = false
	end
end.parse!

def log(msg)
	if $options[:verbose]
		puts msg
	end
end

def err(msg)
	puts msg

	if $options[:fatal]
		exit
	end
end

rename_db = RenameDb.new "protein_rename"
protein_db = ProteinDb.new "human_protein2010.fasta"
xme_list = XMEFilter.new(protein_db, rename_db, ARGV[0])

xme_list.save_database($options[:output])

