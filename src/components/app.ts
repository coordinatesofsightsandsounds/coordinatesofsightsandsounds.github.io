import { Selection, BaseType } from 'd3-selection'
import map from './map'
import video from './video'
import { Emitter } from '../core/eventemitter'


export type Coordinate = {
  degreeLat: number,
  minuteLat: number,
  secondLat: number,
  degreeLng: number,
  minuteLng: number,
  secondLng: number,
}

export type Video = {
  id: number,
  title: string,
  date: string,
  coordinate: Coordinate,
}


export default (app: Selection<BaseType, {}, HTMLElement, any>) => {
  const emitter = new Emitter()

  const videos: { [id: number]: Video } = {
    335598971: {
      id: 335598971,
      title: 'Times Square, New York, NY, USA',
      date: 'December 10, 2018',
      coordinate: {
        degreeLat: 40,
        minuteLat: 45,
        secondLat: 32.2,
        degreeLng: -73,
        minuteLng: 59,
        secondLng: 6.2,
      }
    },
    325310121: {
      id: 325310121,
      title: 'Prospect Park, Brooklyn, NY, USA, March 9, 2019, birds',
      date: 'March 9, 2019',
      coordinate: {
        degreeLat: 40,
        minuteLat: 39,
        secondLat: 30.6,
        degreeLng: -73,
        minuteLng: 58,
        secondLng: 3.7,
      }
    },
  }
  
  const container = app.append('div')

  container.attr('class', 'container')

  map(container, videos, emitter)

  video(container, videos, emitter)
};
