#!/usr/bin/ruby

require 'csv'
require 'set'

require 'read_file.rb'

class RenameDb
	def initialize(filename)
		@protein_rename = {}
		@filename = filename

		load_database
	end

	def load_database
		renames = 0

		file = ReadFile.new(@filename)
		while file.line do
			if file.line =~ /^rename$/
				file.next

				old_names = []
				while file.line !~ /^as$/
					old_names << file.line
					file.next
				end
				file.next
				new_name = file.line
				file.next

				old_names.each do |old_name| 
					@protein_rename[old_name] = new_name
				end

				renames += 1
			end
			file.next
		end
		log "(loaded #{renames} rename sets)"
	end

	def rename(description)
		if @protein_rename.has_key?(description)
			description = @protein_rename[description]
		end

		return description
	end
end
