import { Selection, BaseType, select } from 'd3-selection'
import * as geo from 'd3-geo'
import land from '../assets/ne_110m_land.json'
import { Emitter } from '../core/eventemitter'
import { Video } from './app'
import { degreeToDecimalXY } from '../core/utils'
import { forceSimulation, forceX, forceY, forceCollide } from 'd3-force';
// import * as topojson from 'topojson-client'
// import world from './assets/world-110m.json'
// import eye from './assets/eye.png'

const projection = geo.geoEquirectangular()
// const projection = geo.geoInterruptedHomolosine()
// geo.geoInterruptedSinusoidal();

const markerWidth = 3
const markerHoverWidth = 5


// const simulate = (videos: Video[]) => {
//   const simulation = forceSimulation<Video>(videos)
//     .force('x', forceX<Video>((d) => {
//       return projection(degreeToDecimalXY(d.coordinate))[0]
//     }))
//     .force('y', forceY<Video>((d) => {
//       return projection(degreeToDecimalXY(d.coordinate))[1]
//     }))
//     .force('collide', forceCollide(markerWidth))
//     // .force('collide', (d) => {
//     //   return selection === d.id ? markerHoverWidth : markerWidth
//     // })

//   for (var i = 0; i < 120; ++i) {
//     simulation.tick()
//   }

//   return videos
// }

export default (
  container: Selection<BaseType, {}, HTMLElement, any>,
  videos: Video[],
  bus: Emitter<[number]>
) => {
  
  const simulation = forceSimulation<Video>(videos)
    .force('x', forceX<Video>((d) => {
      return projection(degreeToDecimalXY(d.coordinate))[0]
    }))
    .force('y', forceY<Video>((d) => {
      return projection(degreeToDecimalXY(d.coordinate))[1]
    }))
    .force('collide', forceCollide(markerWidth))

  for (var i = 0; i < 120; ++i) {
    simulation.tick()
  }

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

  svg.append('path')
    .datum({type: 'FeatureCollection', features: land.features})
    .attr('d', geo.geoPath().projection(projection) as any)
    // .attr('d', geo.geoPath().projection(geo.geoMercator()) as any)
    // .attr('d', geo.geoPath().projection(geo.geoNaturalEarth1()) as any)
    // .attr('d', geo.geoPath().projection(geo.geoEqualEarth()) as any)
    .style('stroke-width', '0.6px')


  svg.append('g')
    .selectAll('circle')
    .data(videos)
    .join('circle')
      .attr('cx', d => {
        return d.x
      })
      .attr('cy', d => {
        return d.y
      })
      .attr('r', markerWidth)
      .attr('stroke-width', 1)
      .attr('stroke', 'white')
      .attr('fill', (d) => d.id === -1 ? '#888' : 'red')
      .attr('class', (d) => d.id === -1 ? 'disabled' : '')
      .on('mouseover', function (d) {
        if (d.id !== -1) {
          select(this).attr('fill', '#ad0000')
          // .attr('stroke', '#ffb1b1')
        }
        // select(this).attr('r', markerHoverWidth).attr('stroke-width', 2)
      })
      .on('mouseout', function (d) {
        if (d.id !== -1) {
          select(this).attr('fill', 'red')
          // .attr('stroke', 'white')
        }
        // select(this).attr('r', markerWidth).attr('stroke-width', 1)
      })
      .on('click', ({ id }) => {
        svg.style('display', 'none')
        bus.emit('show-video', id)
      })

  bus.on('hide-video', () => {
    svg.style('display', 'block')
  })
}
