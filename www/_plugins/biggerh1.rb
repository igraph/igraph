# Put the H1 into a div, for bootstrap
module Jekyll
  module BiggerH1Filter
    def biggerh1(input)
      input.gsub(/(<h1[^>]*>.*<\/h1>)/, "<div class=\"page-header\">\\1</div>")
    end
  end
end

Liquid::Template.register_filter(Jekyll::BiggerH1Filter)
