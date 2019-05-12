import { Coordinate } from "../components/app";

export const degreeToDecimalXY = ({
  degreeLat,
  minuteLat,
  secondLat,
  degreeLng,
  minuteLng,
  secondLng,
}: Coordinate): [number, number] => [
  degreeLng + minuteLng/60 + secondLng/3600,
  degreeLat + minuteLat/60 + secondLat/3600,
]

export const degreeToString = ({
  degreeLat,
  minuteLat,
  secondLat,
  degreeLng,
  minuteLng,
  secondLng,
}: Coordinate) => `${Math.abs(degreeLat)}°${minuteLat}'${secondLat}"${degreeLat > 0 ? 'N' : 'S'} ${Math.abs(degreeLng)}°${minuteLng}'${secondLng}"${degreeLng > 0 ? 'E' : 'W'}`
