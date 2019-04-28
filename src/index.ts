import * as geo from 'd3-geo'
import { select } from 'd3-selection'
import * as topojson from 'topojson-client'
// import world from './assets/world-110m.json'
import land from './assets/ne_110m_land.json'
// import eye from './assets/eye.png';



const projection = geo.geoEquirectangular();

// 40°39’30.6”N 73°58’03.7”W
const videos = [
  { lat: 40, lng: -73 }
];

const app = select('#app')

app.append('h1')
  .attr('class', 'title')
  .text('Coordinates of Sights and Sounds')

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
  .attr("d", geo.geoPath().projection(projection) as any)
  // .attr("d", geo.geoPath().projection(geo.geoMercator()) as any)
  // .attr("d", geo.geoPath().projection(geo.geoNaturalEarth1()) as any)
  // .attr("d", geo.geoPath().projection(geo.geoEqualEarth()) as any)
  .style('stroke-width', '0.6px')


svg.append("g")
  .selectAll("circle")
  .data(videos)
  .join("circle")
    .attr("cx", d => projection([d.lng, d.lat])[0])
    .attr("cy", d => projection([d.lng, d.lat])[1])
    .attr("r", 3)

// const markerWidth = 20
// const markerHeight = 20
// svg.append("g")
//   .selectAll("image")
//   .data(videos)
//   .join("image")
//     .attr("xlink:href", eye)
//     .attr("width", `${markerWidth}px`)
//     .attr("height", `${markerHeight}px`)
//     .attr("transform", d => {
//       const [x, y] = projection([d.lng, d.lat])
//       return `translate(${x - (markerHeight / 2)}, ${y - (markerWidth / 2)})`
//     })
