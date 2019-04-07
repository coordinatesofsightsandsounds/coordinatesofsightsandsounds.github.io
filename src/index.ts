import * as geo from 'd3-geo'
import { select } from 'd3-selection'
import * as topojson from 'topojson-client'
// import world from './assets/world-110m.json'
import land from './assets/ne_110m_land.json'


// const svg = document.getElementById('map')
const app = select('#app')

app.append('h1')
  .attr('class', 'title')
  .text('Coordinates of Sites and Sounds')

const svg = app
  .append('div')
  .attr('id', 'map')
  .style('width', 400)
  .append('svg')

svg.attr('viewBox', '0 0 960 600')
  .style('width', '100%')
  .style('height', 'auto')

svg.append("path")
  .datum({type: "FeatureCollection", features: land.features})
  .attr("d", geo.geoPath().projection(geo.geoEquirectangular()) as any)
  // .attr("d", geo.geoPath().projection(geo.geoMercator()) as any)
  // .attr("d", geo.geoPath().projection(geo.geoNaturalEarth1()) as any)
  // .attr("d", geo.geoPath().projection(geo.geoEqualEarth()) as any)
  .style('fill', 'white')
  .style('stroke', '#000')
  .style('stroke-width', '0.5px');

// svg.selectAll("path")
//   .data(land)
//   .enter().append("path")
//     .attr("d", geo.geoPath());

// svg.append('path')
//   .datum(land)
//   .enter().append('path')
//     .attr('d', geo.geoPath().projection(geo.geoMercator()))

// svg.selectAll("path")
//   // .data(topojson.feature(world, world.objects.land).features)
//   // .datum(topojson.mesh(world, world.objects.land))
//   .enter().append("path")
//   .attr("d", path);
