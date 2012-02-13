#!/usr/bin/ruby

require 'csv'
require 'set'

require 'read_file.rb'

class Fraction
	def initialize(filename)
		@protein_fraction = {}
		@filename = filename

		load_database
	end

	def load_database
		log "loading fraction list #{@filename} ..."

		CSV.open(@filename, 'r') do |row|
			desc = row[1].strip
			fraction = row[4].strip

			@protein_fraction[desc] = fraction
		end
	end

	def [](key)
		@protein_fraction[key]
	end
end
