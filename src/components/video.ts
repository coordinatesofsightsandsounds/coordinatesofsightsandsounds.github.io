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
      .text('×')
      // .text('x')
      .on('click', () => {
        bus.emit('hide-video', id)
        video.remove()
      })

    // video.append('h3').text(videos[id].title)

    video.append('h3').text(videos[id].date).attr('class', 'date')
    
    video.append('div').attr('id', 'player')

    // video.append('h3').text('So and So').attr('class', 'attribution')


    // TODO - might be possible to initialize the player and change the id
    var player = new Player('player', {
      id,
      width: 640,
      loop: true
    })

    player.play();
  })
}
