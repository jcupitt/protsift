#!/usr/bin/ruby

require 'csv'
require 'set'

#glutathione S-transferase alpha
#
#>gi|4504173|ref|NP_001503.1| glutathione S-transferase alpha 4 [Homo sapiens]
#>gi|215276987|ref|NP_000837.3| glutathione S-transferase alpha 2 [Homo sapiens]
#>gi|24430144|ref|NP_000838.3| glutathione S-transferase alpha 3 [Homo sapiens]
#>gi|24308514|ref|NP_714543.1| glutathione S-transferase alpha 5 [Homo sapiens]
#>gi|22091454|ref|NP_665683.1| glutathione S-transferase alpha 1 [Homo sapiens]

class XMEFilter
	def initialize(protein_db, rename_db, filename)
		@protein_db = protein_db
		@rename_db = rename_db
		@filename = filename

		# hash from description to a set of NP numbers
		@family_table = Hash.new {|h, k| h[k] = Set.new}

		# hash from a description to the description split into a set of
		# search terms
		@search_terms = {}

		load_database
	end

	def [](key)
		@family_table[key]
	end

	def each_desc
		@family_table.each_key do |desc|
			yield desc
		end
	end

	# does one of our short-form renamed descriptions (eg. "cytochrome P450
	# 2A6/2A7/2A13") match a long-form description from the fasta db (eg.
	# "cytochrome P450 2A7 isoform 2")
	def is_xme(desc, description)
		terms = @search_terms[desc]

		match = true
		terms.each do |term|
			word_match = false
			term.each do |word|
				if description.include?(word)
					word_match = true
					break
				end
			end

			if not word_match
				match = false
				break
			end
		end

		return match
	end

	# given a long-form protein description, find the corresponding xme desc, 
	# or nil if none
	def find_desc(description)
		description = @rename_db.rename(description)

		each_desc do |desc|
			if is_xme(desc, description)
				return desc
			end
		end

		return nil
	end

	def load_database
		log "loading XME database #{@filename} ..."

		linen = 0
		CSV.open(@filename, 'r') do |row|
			# all useful ones should begin with "NP", or have an empty column 1
			# and something in column 2
			nprow = row[0] =~ /^NP/
			descrow = (row[0] == nil and row[1] != nil)
			if not nprow and not descrow
				next
			end
	
			nps = row[0]
			desc = row[1]
			if nps == nil
				nps = ""
			end

			# stray spaces appear in the nps
			nps = nps.strip
			nps = nps.split(",")
			nps = nps.map {|x| x.strip}

			# the desc line often has stray spaces at the end
			desc = desc.strip

			@family_table[desc] += nps
		end

		# make all the search queries
		@family_table.each do |desc, nps|
			terms = desc.split(" ")
			terms.map! {|term| term.split("/")}

			@search_terms[desc] = terms
		end

		@family_table.each do |desc, nps|
			found = false
			@protein_db.each do |protein|
				description = @rename_db.rename(protein.description)
				if is_xme(desc, description)
					nps << protein.ncbi
					found = true
				end
			end
			if not found
				err "no protein found for #{desc}"
			end
		end
	end

	def save_database(filename)
		log "writing #{filename} ..."
		CSV.open(filename, 'w') do |write|
			row = []
			row << "ncbi"
			row << "description"
			write << row

			@family_table.keys.sort.each do |desc|
				nps = @family_table[desc]

				row = []
				row << nps.to_a.join(", ")
				row << desc
				write << row
			end
		end
	end
end

