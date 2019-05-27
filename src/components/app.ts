import { Selection, BaseType } from 'd3-selection'
import map from './map'
import video from './video'
import { Emitter } from '../core/eventemitter'
import { decimalToDegree } from '../core/utils';
import { SimulationNodeDatum } from 'd3-force';


export type Coordinate = {
  degreeLat: number,
  minuteLat: number,
  secondLat: number,
  degreeLng: number,
  minuteLng: number,
  secondLng: number,
}

export type Video = SimulationNodeDatum & {
  id: number,
  title: string,
  date: string,
  coordinate: Coordinate,
}


export default (app: Selection<BaseType, {}, HTMLElement, any>) => {
  const emitter = new Emitter()

  const videos: { [id: number]: Video } = {
    335594697: {
      id: 335594697,
      title: 'Coney Island, Brooklyn, NY, USA, June 16, 2018, under jetty (light organ)',
      date: 'June 16, 2018',
      coordinate: decimalToDegree([-73.983913, 40.572231])
    },
    335595236: {
      id: 335595236,
      title: 'Lake Sacandaga outlet, Speculator, NY, USA, July 7, 2018, moon ',
      date: 'June 16, 2018',
      coordinate: decimalToDegree([-74.404956, 43.474052])
    },
    335596045: {
      id: 335596045,
      title: 'River Road, NJ, USA, July 9, 2018, fireflies',
      date: 'July 9, 2018',
      coordinate: decimalToDegree([-74.857114, 40.283332])
    },
    335597358: {
      id: 335597358,
      title: 'Collioure, France',
      date: 'August 12, 2018',
      coordinate: decimalToDegree([3.086859, 42.527878])
    },
    335598158: {
      id: 335598158,
      title: 'Prospect Park, Brooklyn, NY, USA, November 10, 2018, tunnel',
      date: 'November 10, 2018',
      coordinate: decimalToDegree([-73.964821, 40.659920])
    },
    335598403: {
      id: 335598403,
      title: 'Roosevelt Island Tramway, New York, NY, USA, November 15, 2018',
      date: 'November 15, 2018',
      coordinate: decimalToDegree([-73.957244, 40.758128])
    },
    335598971: {
      id: 335598971,
      title: 'Times Square, New York, NY, USA, December 10, 2018',
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
    10: {
      id: -1,
      title: 'Fort Greene Park, Brooklyn, NY, USA, March 15, 2019, top of the hill',
      date: 'March 15, 2019',
      coordinate: {
        degreeLat: 40,
        minuteLat: 41,
        secondLat: 30.8,
        degreeLng: -73,
        minuteLng: 58,
        secondLng: 32.3,
      }
    },
    11: {
      id: -1,
      title: 'East River Esplanade, 62nd Street, New York, NY, USA, March 21, 2019',
      date: 'March 21, 2019',
      coordinate: {
        degreeLat: 40,
        minuteLat: 45,
        secondLat: 36.1,
        degreeLng: -73,
        minuteLng: 57,
        secondLng: 25.6,
      }
    },
    12: {
      id: -1,
      title: 'Guggenheim Museum,  New York, NY, USA, March 23, 2019',
      date: 'March 23, 2019',
      coordinate: {
        degreeLat: 40,
        minuteLat: 46,
        secondLat: 59.5,
        degreeLng: -73,
        minuteLng: 57,
        secondLng: 32.4,
      }
    },
    13: {
      id: -1,
      title: 'Bucks County Playhouse, New Hope, NJ, USA',
      date: 'April 21, 2019',
      coordinate: decimalToDegree([-74.950765, 40.362344])
    },
    14: {
      id: -1,
      title: 'Delaware and Raritan Canal State Park Trail, Lambertville, NJ, USA, April 21, 2019',
      date: 'April 21, 2019',
      coordinate: decimalToDegree([-74.945191, 40.363840])
    },
    15: {
      id: -1,
      title: 'Chaim Baier Music Island, Prospect Park, Brooklyn, NY, USA, May 4, 2019',
      date: 'May 4, 2019',
      coordinate: decimalToDegree([-73.965041, 40.658142])
    },
  }
  
  const container = app.append('div')

  container.attr('class', 'container')

  map(container, Object.keys(videos).map((id) => videos[id]), emitter)

  video(container, videos, emitter)
};
