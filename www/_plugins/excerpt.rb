# Include everything from the title until the first header in the excerpt
# http://tech.pro/tutorial/1299/getting-started-with-jekyll-plugins
module Jekyll
  module ExcerptFilter
    def extract_excerpt(input, url)
      input=input.split("</h1>", 2)[1]
      re = %r{<!--more(?<read_more_message> ...*)?-->}
      m = input.match re

      if m.nil?
        input
      elsif m[:read_more_message].nil?
        %|#{input.split(m[0]).first.strip}<p><a href="#{url}">More &#8594;</a></p>|
      else
        %|#{input.split(m[0]).first.strip}<p><a href="#{url}">#{m[:read_more_message].strip}</a></p>|
      end
    end
  end
end

Liquid::Template.register_filter(Jekyll::ExcerptFilter)
