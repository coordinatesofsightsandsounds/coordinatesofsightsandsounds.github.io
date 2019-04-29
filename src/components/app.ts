import { Selection, BaseType } from 'd3-selection';
import map from './map';


export default (app: Selection<BaseType, {}, HTMLElement, any>) => {
  app.append('h1')
    .attr('class', 'title')
    .text('Coordinates of Sights and Sounds')

  map(app);
};
