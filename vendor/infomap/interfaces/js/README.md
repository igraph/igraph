# Infomap

This is Infomap compiled as a web worker with Emscripten.

Infomap is a network clustering algorithm based on the [Map equation](//www.mapequation.org/publications.html#Rosvall-Axelsson-Bergstrom-2009-Map-equation).

This package is used in [Infomap Online](//www.mapequation.org/infomap/).

## Installing

To install, run

```shell
npm install @mapequation/infomap
```

## Usage

If you use ES modules, import the package like this

```javascript
import Infomap from "@mapequation/infomap";

let network = `#source target [weight]
0 1
0 2
0 3
1 0
1 2
2 1
2 0
3 0
3 4
3 5
4 3
4 5
5 4
5 3`;

let infomap = new Infomap()
  .on("data", (data) => console.log(data))
  .on("error", (err) => console.warn(err))
  .on("finished", (data) => console.log(data));

infomap.run(network, "--two-level");
```

If you use a CDN, for example JSDelivr, `Infomap` is exported as `window.infomap.default`.

For example:

```html
<!doctype html>
<html>
    <head>
        <script type="text/javascript" src="//cdn.jsdelivr.net/npm/@mapequation/infomap@latest/dist/index.min.js"></script>
    </head>
    <body>
        <script type="text/javascript">
            const Infomap = window.infomap.default;

            let network = "#--- as above! ---";

            let infomap = new Infomap()
                .on("data", data => console.log(data))
                .on("error", err => console.warn(err))
                .on("finished", data => console.log(data));

            infomap.run(network, "--two-level");
        </script>
</html>
```

## Feedback

If you have any questions, suggestions or issues regarding the software, please add them to [GitHub issues](//github.com/mapequation/infomap/issues).

## Authors

Daniel Edler, Anton Holmgren, Martin Rosvall

For contact information, see [mapequation.org/about.html](//www.mapequation.org/about.html).

## Terms of use

Infomap is released under a dual licence.

To give everyone maximum freedom to make use of Infomap and derivative works,
we make the code open source under the GNU Affero General Public License version 3
or any later version (see LICENSE_AGPLv3.txt).

For a non-copyleft license, please contact us.
