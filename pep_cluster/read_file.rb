#!/usr/bin/ruby

# read a file a line at a time, with a lookahead token
# handy for parsing!
class ReadFile 
	def initialize(filename)
		@filename = filename
		@file = File.new(filename)
		@line = @file.gets
	end

	def line
		if @line == nil
			nil
		else
			@line.chomp
		end
	end

	def next
		@line = @file.gets
	end
end

