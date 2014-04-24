# Rewrite old docs links with an absolute url back to SF
module Jekyll
  module RedirDocsFilter
    def redir_docs(input)
      input.gsub("href=\"doc-", "href=\"http://igraph.sourceforge.net/doc-")
    end
  end
end

Liquid::Template.register_filter(Jekyll::RedirDocsFilter)
