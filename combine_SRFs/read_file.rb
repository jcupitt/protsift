#!/usr/bin/ruby

# read a file a line at a time, with a lookahead token
# handy for parsing!

class ReadFile 
	def initialize(filename)
		@filename = filename

		log "loading #{@filename} ..."

		@file = File.new(filename)

		self.next
	end

	def line
		if @line == nil
			nil
		else
			# remove trailing whitespace too, very annoying and confusing in
			# config files
			@line.chomp.rstrip
		end
	end

	def next
        # skip lines starting with a hash
        begin
		    @line = @file.gets
        while @line[0] == "#"
	end
end

