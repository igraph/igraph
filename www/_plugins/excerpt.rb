# Include everything from the title until the first header in the excerpt
module Jekyll
  module ExcerptFilter
    def extract_excerpt(input)
      input.split("<h2")[0].split("</h1>")[1]
    end
  end
end
 
Liquid::Template.register_filter(Jekyll::ExcerptFilter)
