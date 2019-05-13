import { Selection, BaseType, select } from 'd3-selection';
import * as geo from 'd3-geo'
import land from '../assets/ne_110m_land.json'
import { Emitter } from '../core/eventemitter';
import { Video } from './app';
import { degreeToDecimalXY } from '../core/utils';
// import * as topojson from 'topojson-client'
// import world from './assets/world-110m.json'
// import eye from './assets/eye.png';

const projection = geo.geoEquirectangular();

const markerWidth = 3
const markerHoverWidth = 5

export default (
  container: Selection<BaseType, {}, HTMLElement, any>,
  videos: { [id: number]: Video },
  bus: Emitter<[number]>
) => {
  container.append('h1')
    .attr('class', 'title')
    .text('Coordinates of Sights and Sounds')

  const svg = container
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
    .data(Object.keys(videos).map((id) => videos[id]))
    .join("circle")
      .attr("cx", d => projection(degreeToDecimalXY(d.coordinate))[0])
      .attr("cy", d => projection(degreeToDecimalXY(d.coordinate))[1])
      .attr("r", markerWidth)
      .attr("stroke-width", 1)
      .on("mouseover", function () {
        select(this).attr("r", markerHoverWidth).attr("stroke-width", 2)
      })
      .on("mouseout", function () {
        select(this).attr("r", markerWidth).attr("stroke-width", 1)
      })
      .on("click", ({ id }) => {
        svg.style("display", "none")
        bus.emit("show-video", id)
      })

  bus.on('hide-video', () => {
    svg.style("display", "block")
  })

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
}
