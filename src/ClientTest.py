import websockets
import asyncio

URL = "ws://127.0.0.1:8080/websocket/"
NAME = "wyz"


async def test(url):
    async with websockets.connect(url) as websocket:
        await websocket.send("hello")
        while True:
            recv = await websocket.recv()
            print(recv)


if __name__ == '__main__':
    asyncio.get_event_loop().run_until_complete(test(URL + NAME))
