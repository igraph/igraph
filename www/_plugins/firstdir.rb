# Extract the top directory from a URL
module Jekyll
  module FirstdirFilter
    def firstdir(input)
      input.gsub(/(\/[^\/]*\/?).*$/, '\1')
    end
  end
end

Liquid::Template.register_filter(Jekyll::FirstdirFilter)
