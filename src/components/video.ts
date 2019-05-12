import { Selection, BaseType } from 'd3-selection';
import Player from '@vimeo/player';
import { Emitter } from '../core/eventemitter';
import { Video } from './app';
import { degreeToString } from '../core/utils';

export default (
  container: Selection<BaseType, {}, HTMLElement, any>,
  videos: { [id: number]: Video },
  bus: Emitter<[number]>
) => {
  bus.on('show-video', (id) => {
    const video = container
      .append('div')
      .attr('id', 'video')

    video.append('h2').text(degreeToString(videos[id].coordinate))
      .append('button')
      .attr('class', 'close')
      .text('X')
      .on('click', () => video.remove())

    video.append('h3').text(videos[id].title)

    video.append('h3').text(videos[id].date)
    
    video.append('div').attr('id', 'player')


    // TODO - might be possible to initialize the player and change the id
    var player = new Player('player', {
      id,
      width: 640,
      loop: true
    })
  
    player.setVolume(0);
  
    player.on('player', function() {
      console.log('played the video!')
    })
  })
}
